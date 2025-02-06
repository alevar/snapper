import re
import sys
import pysam
import argparse

from .utils.common import *
from .classes.txgroup import Transcriptome

class Snapper:
    def __init__(self, args):
        # input/output files
        self.sam = args.sam
        self.reference = args.reference
        self.output = args.output
        
        # scoring parameters
        self.qry_intron_match_score = args.qry_intron_match_score
        self.trg_pos_match_score = args.trg_pos_match_score
        self.trg_pos_mismatch_score = args.trg_pos_mismatch_score

    def run(self):
        ref_tome = Transcriptome()
        ref_tome.build_from_file(self.reference)
        self.adjust_cigar_with_pysam(self.sam, self.output, ref_tome)

    def adjust_cigar_with_pysam(self, input_bam, output_bam, ref_tome):
        """Precisely align reference intron boundaries while preserving other features"""
        
        def parse_cigar(cigar_string):
            """Parse CIGAR string into operations and lengths"""
            return [(int(length), op) for length, op in re.findall(r'(\d+)(\D)', cigar_string)]
        
        def cigar_to_blocks(cigar_ops, query_pos=0, target_pos=0):
            """Convert CIGAR operations to blocks of matching/mismatching regions"""
            blocks = []
            current_pos = {'query': query_pos, 'target': target_pos}
            
            for length, op in cigar_ops:
                if op == 'N':
                    blocks.append(('N', length, current_pos['query'], current_pos['target']))
                    current_pos['target'] += length
                else:
                    if op in ('M', 'X', '='):
                        blocks.append((op, length, current_pos['query'], current_pos['target']))
                        current_pos['query'] += length
                        current_pos['target'] += length
                    elif op == 'I':
                        blocks.append((op, length, current_pos['query'], current_pos['target']))
                        current_pos['query'] += length
                    elif op == 'D':
                        blocks.append((op, length, current_pos['query'], current_pos['target']))
                        current_pos['target'] += length
                    elif op == 'S':
                        blocks.append((op, length, current_pos['query'], current_pos['target']))
                        current_pos['query'] += length
            return blocks

        def score_alignment(blocks, intron_positions, qry_match_score, pos_match_score, mismatch_score):
            """Score alignment based on intron positions and matching positions"""
            score = 0
            for block in blocks:
                op, length, qpos, tpos = block
                if op == 'N':
                    if qpos in intron_positions:
                        score += qry_match_score * length
                elif op in ('M', 'X', '='):
                    for i in range(length):
                        if (qpos + i) == (tpos + i):
                            score += pos_match_score
                        else:
                            score += mismatch_score
            return score

        def blocks_to_cigar(blocks):
            """Convert blocks back to CIGAR string"""
            cigar_ops = []
            current_op = None
            current_length = 0
            
            for block in blocks:
                op, length = block[0], block[1]
                if op == current_op:
                    current_length += length
                else:
                    if current_op:
                        cigar_ops.append(f"{current_length}{current_op}")
                    current_op = op
                    current_length = length
            
            if current_op:
                cigar_ops.append(f"{current_length}{current_op}")
            
            return "".join(cigar_ops)

        def get_intron_positions(exons):
            """Get intron positions from exon coordinates"""
            if len(exons) <= 1:
                return []
            positions = []
            pos = 0
            for start, end in exons[:-1]:
                pos += (end - start)
                positions.append(pos - 1)  # -1 because we want 0-based position
            return positions

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
                    exons = [(x[0]-1, x[1]-1) for x in transcript.get_exons()]
                    intron_positions = get_intron_positions(exons)
                    
                    # Parse CIGAR and convert to blocks
                    cigar_ops = parse_cigar(read.cigarstring)
                    blocks = cigar_to_blocks(cigar_ops)
                    
                    # Score the alignment and adjust blocks if needed
                    best_score = score_alignment(blocks, intron_positions, 
                                              self.qry_intron_match_score,
                                              self.trg_pos_match_score, 
                                              self.trg_pos_mismatch_score)
                    best_blocks = blocks
                    
                    # Convert back to CIGAR string
                    new_cigar = blocks_to_cigar(best_blocks)
                    
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