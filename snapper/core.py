import os
import re
import sys
import copy
import pysam
import shutil
import argparse
import subprocess

from itertools import combinations

from .utils.common import *

from .classes.txgroup import Transcriptome, Gene, Bundle
from .classes.transcript import Transcript, Object

class Snapper:
    def __init__(self, args):         
        # INPUT FILES
        self.sam = args.sam
        self.reference = args.reference
        self.output = args.output

        self.qry_intron_match_score = args.qry_intron_match_score
        self.trg_pos_match_score = args.trg_pos_match_score
        self.trg_pos_mismatch_score = args.trg_pos_mismatch_score

    def run(self):
        # start by building transcriptomes for reference and target
        ref_tome = Transcriptome()
        ref_tome.build_from_file(self.reference)
        
        self.adjust_cigar_with_pysam(self.sam, self.output, ref_tome)

    def adjust_cigar_with_pysam(self, input_bam, output_bam, ref_tome):
        """Precisely align reference intron boundaries while preserving other features"""
        
        def parse_cigar(cigar):
            return [(int(length), op) for length, op in re.findall(r'(\d+)(\D)', cigar)]

        def build_cigar(ops):
            """Convert list of CIGAR tuples to CIGAR string"""
            return ''.join(f"{length}{op}" for op, length in ops)
        
        def explode_cigar(ops):
            """
            Explode CIGAR operations into single-base pair elements.
            
            Args:
                ops: List of CIGAR tuples in format [(length, operation)]
                
            Returns:
                List of single-base pair operations
                
            Example:
                >>> explode_cigar([(3, 'M'), (2, 'I')])
                ['M','M','M','I', 'I']
            """
            exploded = []
            for length, op in ops:
                exploded.extend([op] * length)
            return exploded
        
        def exons_to_introns(exons):
            ref_intron_vec = []
            if len(exons) > 1:
                pos = 0
                for exon in exons:
                    pos += (exon[1] - exon[0])
                    ref_intron_vec.append(pos)
                ref_intron_vec.pop() # remove last position since it's not an intron
            return ref_intron_vec

        def extract_cigar_data(cigar, qry_intron_positions):
            # Separate cigar into blocks of N and everything else
            aln_ops = []  # Stores alignment operations
            intron_blocks = []  # Stores intronic blocks

            qry_pos = 0  # Tracks the position in the query sequence - query position will never change since not affected by introns
            trg_pos = 0  # Tracks the position in the reference sequence - will depend on the intron configuration

            # Temporary storage for the current block
            current_block = []
            is_intron = False
            
            intronic_positions = dict() # positions in the original cigar that are followed by introns

            for op in cigar:
                if op == 'N':
                    # If we encounter 'N', it's an intron block
                    if not is_intron and current_block:
                        # If we were in an alignment block, store it
                        aln_ops.extend(current_block)
                        current_block = []
                    if not is_intron:
                        ipos = len(aln_ops)
                        intronic_positions.setdefault(ipos, 0) # set to 0 to indicate that this is a intorn on target, but if already set to 1, then it's a query intron and stays that way
                    is_intron = True
                    current_block.append('N')
                    trg_pos += 1  # 'N' advances the reference position but not the query position
                else:
                    if qry_pos in qry_intron_positions:
                        ipos = len(aln_ops)
                        if not is_intron:
                            ipos += len(current_block)
                        intronic_positions[ipos] = 1 # set to 1 to indicate that this is a query position
                    # If we encounter anything else, it's an alignment operation
                    if is_intron and current_block:
                        # If we were in an intron block, store it
                        intron_blocks.append(current_block)
                        current_block = []
                    is_intron = False
                    current_block.append((op, qry_pos, trg_pos))
                    if op in ['M', 'X', '=']:  # These operations advance both query and reference positions
                        qry_pos += 1
                        trg_pos += 1
                    elif op == 'S': # Soft clipping advances query position only
                        qry_pos += 1
                    elif op == 'I':  # Insertion advances query position only
                        qry_pos += 1
                    elif op == 'D':  # Deletion advances reference position only
                        trg_pos += 1

            # Append the last block if it exists
            if current_block:
                if is_intron:
                    intron_blocks.append(current_block)
                else:
                    aln_ops.extend(current_block)
                    
            return aln_ops, intron_blocks, intronic_positions

        def generate_cigar_variations(ops, intronic_positions, intron_blocks):
            """
            Generate all possible cigar variations by inserting introns at specified positions.

            Args:
                ops (list): The input list of operations.
                intronic_positions (list): List of positions in the list where insertion can happen.
                intron_blocks (list): List of strings to insert in the order specified.

            Returns:
                list: A list of all possible variations of the list `s`.
            """
            if len(intron_blocks) > len(intronic_positions):
                raise ValueError("Number of intron_blocks cannot be greater than the number of intronic positions.")

            cigars = []

            # Generate all combinations of positions to insert the strings
            for pos_comb in combinations(intronic_positions, len(intron_blocks)):
                temp_s = ops[:]
                offset = 0  # To adjust the insertion index as the list changes
                for iblock, pos in zip(intron_blocks, pos_comb):
                    temp_s = temp_s[:pos + offset] + [(x,intronic_positions[pos]) for x in iblock] + temp_s[pos + offset:]
                    offset += len(iblock)  # Update offset after insertion
                cigars.append(temp_s)

            return cigars
        
        def score_cigar(cigar, qry_intron_match_score, trg_pos_match_score, trg_pos_mismatch_score):
            score = 0

            qry_pos = 0  # Tracks the position in the query sequence - query position will never change since not affected by introns
            trg_pos = 0  # Tracks the position in the reference sequence - will depend on the intron configuration

            is_intron = False

            for op in cigar:
                if op[0] == 'N':
                    if not is_intron: # first intronic position - check the score
                        score += op[1] * qry_intron_match_score # op[1] set to 1 if matches the query intron position, 0 otherwise
                    is_intron = True
                    trg_pos += 1  # 'N' advances the reference position but not the query position
                else:
                    is_intron = False
                    if trg_pos == op[2]: # matching position
                        score += trg_pos_match_score
                    else:
                        score += trg_pos_mismatch_score
                    if op[0] in ['M', 'X', '=']:  # These operations advance both query and reference positions
                        qry_pos += 1
                        trg_pos += 1
                    elif op[0] == 'S': # Soft clipping advances query position only
                        qry_pos += 1
                    elif op[0] == 'I':  # Insertion advances query position only
                        qry_pos += 1
                    elif op[0] == 'D':  # Deletion advances reference position only
                        trg_pos += 1
                        
            return score
        
        def merge_ops(ops):
            merged_cigar = []
            count = 1
            for i in range(1, len(ops)):
                if ops[i][0] == ops[i-1][0]:
                    count += 1
                else:
                    merged_cigar.append((ops[i-1][0], count))
                    count = 1
            merged_cigar.append((ops[-1][0], count))
            return merged_cigar
        
        def adjust_cigar(original_cigar, exons):
            """Adjust CIGAR to match reference exon boundaries exactly"""

            # Parse and explode original CIGAR
            raw_ops = parse_cigar(original_cigar)
            ops = explode_cigar(raw_ops)

            exons = [(x[0]-1, x[1]-1) for x in exons]
            qry_intron_positions = exons_to_introns(exons)

            aln_ops, intron_blocks, trg_intronic_positions = extract_cigar_data(ops, qry_intron_positions)

            # generate all possible variations of the alignment operations
            cigar_variants = generate_cigar_variations(aln_ops, trg_intronic_positions, intron_blocks)
            
            best_cigar = max(cigar_variants, key=lambda x: score_cigar(x, self.qry_intron_match_score, self.trg_pos_match_score, self.trg_pos_mismatch_score))

            # Merge operations back to standard CIGAR format
            res_cigar = build_cigar(merge_ops(best_cigar))
            return res_cigar

        # File handling with pysam
        input_mode = 'rb' if input_bam.endswith('.bam') else 'r'
        output_mode = 'wb' if output_bam.endswith('.bam') else 'wh'

        with pysam.AlignmentFile(input_bam, input_mode) as infile, \
            pysam.AlignmentFile(output_bam, output_mode, template=infile) as outfile:

            for read in infile:
                if read.is_unmapped or not read.cigarstring:
                    outfile.write(read)
                    continue

                try:
                    transcript = ref_tome.get_by_tid(read.query_name)
                    
                    new_cigar = adjust_cigar(read.cigarstring, transcript.get_exons())
                    
                    # Create modified read
                    modified_read = pysam.AlignedSegment(outfile.header)
                    modified_read.set_tags(read.get_tags())
                    modified_read.query_name = read.query_name
                    modified_read.query_sequence = read.query_sequence
                    modified_read.flag = read.flag
                    modified_read.reference_id = read.reference_id
                    modified_read.reference_start = read.reference_start
                    modified_read.mapping_quality = read.mapping_quality
                    modified_read.cigarstring = new_cigar
                    modified_read.query_qualities = read.query_qualities
                    modified_read.next_reference_id = read.next_reference_id
                    modified_read.next_reference_start = read.next_reference_start
                    modified_read.template_length = read.template_length

                    outfile.write(modified_read)

                except Exception as e:
                    sys.stderr.write(f"Error processing {read.query_name}: {str(e)}\n")
                    outfile.write(read)
        
        
def main():
    parser = argparse.ArgumentParser(description="Correct Intron shifts in alignments via reference annotation")

    parser.add_argument('-s', '--sam', required=True, type=str, help='Path to the sam alignment file')
    parser.add_argument('-r', '--reference', required=True, type=str, help='Path to the reference annotation')
    parser.add_argument('-o', '--output', type=str, help='Path to the output GTF file')
    
    parser.add_argument('--qry_intron_match_score', type=int, default=10, help='Score for matching query introns')
    parser.add_argument('--trg_pos_match_score', type=int, default=1, help='Score for matching target positions')
    parser.add_argument('--trg_pos_mismatch_score', type=int, default=-1, help='Score for mismatching target positions')

    args = parser.parse_args()

    snapper = Snapper(args)
    snapper.run()

if __name__ == "__main__":
    main()