from intervaltree import Interval, IntervalTree
from typing import Iterator, List, Callable
import copy
import os

from pyfaidx import Fasta

from ..classes.treader import TReader
from ..classes.transcript import Transcript, Object
from ..utils.common import *

class TXGroup:
    """
    Represents a general class of a group of objects.
    In this general implementation any object can be added to the group.

    Attributes:
        objects (list): A list of objects.
        tid_map (dict): A mapping of transcript IDs to positions in the objects list.
    """
    def __init__(self):
        self.objects = list()
        self.tid_map = dict() # transcript_id to position in objects list
        
    def __contains__(self,tid) -> bool:
        return tid in self.tid_map
    
    def __len__(self) -> int:
        return len(self.objects)

    def clear(self):
        """
        Clear the objects and tid_map.
        """
        self.objects = list()
        self.tid_map = dict()
        
    def add_object(self,obj) -> int:
        """
        Add an object to the TXGroup. In this general implementation any object can be added to the group.

        Args:
            obj: The object to add.

        Returns:
            int: The index of the added object in the objects list.
        """
        idx = None
        if obj.get_type() in [Types.Transcript, Types.MRNA, Types.Exon, Types.CDS]:
            idx = self.tid_map.setdefault(obj.get_tid(),len(self.objects))
            if len(self.objects) == idx:
                self.objects.append(Transcript(obj))
            else:
                self.objects[idx].merge(obj)
        elif obj.get_type() in [Types.Bundle, Types.Gene]:
            for tx in obj.transcript_it():
                idx = self.tid_map.setdefault(tx.get_tid(),len(self.objects))
                assert len(self.objects) == idx,"Transcript ID already exists in the objects list"
                self.objects.append(tx)
        elif obj.get_type() in [Types.Other]:
            pass
        else:
            raise Exception("Wrong object type provided to the add_object methods of TXGroup")
        
        return idx
    
    def idx(self,tid) -> int:
        """
        Get the index of an object in the TXGroup by transcript ID.

        Args:
            tid: The transcript ID of the object.

        Returns:
            int: The index of the object in the objects list.
        """
        return self.tid_map.get(tid,None)
    
    def __getitem__(self,idx):
        """
        Get an object from the TXGroup by index.

        Args:
            idx: The index of the object.

        Returns:
            object: The object at the specified index.
        """
        return self.objects[idx]
    
    def get_by_tid(self,tid) -> Object:
        """
        Get an object from the TXGroup by transcript ID.

        Args:
            tid: The transcript ID of the object.

        Returns:
            object: The object with the specified transcript ID.
        """
        idx = self.idx(tid)
        if idx is not None:
            return self.objects[idx]
        return None
        
    def object_it(self, func: Callable=None) -> Iterator[Object]:
        """
        Iterate over all objects in the TXGroup.
        Specific implementations of the TXGroup may provide specialized implementations such as transcdrip_it which yield subsets of the objects.
        This one should not be overwritten.

        Args:
            func (Callable, optional): A function to apply to each object before yielding it. If function returns False, the object is skipped.keep
        """
        if func is None:
            yield from self.objects
        else:
            for obj in self.objects:
                if func(obj):
                    yield obj

    def is_empty(self) -> bool:
        """
        Check if the TXGroup is empty.

        Returns:
            bool: True if the TXGroup is empty, False otherwise.
        """
        return len(self.objects) == 0
    
    def group_by(self,by) -> None:
        """
        Sorts objects internally according to the desired criteria and then yields groups of them

        Args:
            by (str, optional): The criteria to group by..
        """

        self.sort(by)

        prev_key = None
        group = TXGroup()
        no_cds = []
        for obj in self.objects:
            try:
                if not obj.has_cds():
                    no_cds.append(obj)
                    continue
            except:
                pass

            key = obj._getattrs(by)
            if key != prev_key:
                if not group.is_empty():
                    yield prev_key, group
                group = TXGroup()
            group.add_object(obj)
            prev_key = key
        if not group.is_empty():
            yield prev_key, group
        for obj2 in no_cds:
            group = TXGroup()
            group.add_object(obj2)
            yield obj2._getattrs(by), group

    def sort(self,by) -> None:
        """
        Sort the objects in the TXGroup.
        Args:
            by (str or list, required): The criteria to sort by.
        """
        cmp = lambda obj: obj._getattrs(by)
        self._sort(cmp)

    def size(self) -> int:
        return len(self.objects)
    
    def _sort(self,cmp) -> None:
        """
        Sort the objects in the TXGroup using a custom comparison function.

        Args:
            cmp: The comparison function to determine the sorting order.
        """
        self.objects.sort(key=cmp)
        self.reindex()


    def to_gtf(self):
        """
        Convert the TXGroup to GTF format.

        Returns:
            str: The TXGroup objects formatted as a GTF string.
        """
        res = ""
        for obj in self.objects:
            res+=obj.to_gtf()+"\n"
        res.rstrip("\n")
        return res
    
    def to_gff(self):
        """
        Convert the TXGroup to GFF format.

        Returns:
            str: The TXGroup objects formatted as a GFF string.
        """
        res = ""
        for obj in self.objects:
            res+=obj.to_gtf()+"\n"
        res.rstrip("\n")
        return res

    # for every object in the object list reconstruct the tidmap
    def reindex(self):
        """
        Reconstruct the tid_map based on the objects in the TXGroup.
        """
        self.tid_map.clear()
        idx = 0
        for obj in self.objects:
            assert obj.get_tid() not in self.tid_map,"duplicate object IDs: "+obj.get_tid()
            self.tid_map[obj.get_tid()] = idx
            idx+=1
    
    __repr__ = to_gtf

    __str__ = __repr__

class Bundle (TXGroup,Object):
    """
    Specialization designed to describe an Object which is also a Group of Objects 
    such as OverlapBundle (transcripts with overlapping coordinates)
    or Genes (transcripts with the same geneID). 
    Inherits also from the common object traits such as start,end,overlaps, etc

    Unlike TXGroup, since Bundle 9inherits from Object, it has start and end and is guaranteed to exist on the same seqid and strand
    """
    def __init__(self):
        TXGroup.__init__(self)
        Object.__init__(self)

        self.intervals = IntervalTree() # union of all exons in the locus (minus the introns)

    def  add_object(self,obj: Object) -> bool:
        """
        Add an object to the Bundle.
        Only adds if satisfies basic Bundle constraints:
        1. all objects share seqid
        2. all objects share strand
        Updates start/end

        Args:
            obj (Object): The object to add.

        Returns:
            bool: True if the object was successfully added, False otherwise.
        """
        if self.seqid is None:
            self.seqid = obj.seqid
        if self.strand is None:
            self.strand = obj.strand

        if self.seqid != obj.get_seqid() or self.strand != obj.get_strand():
            return False

        idx = TXGroup.add_object(self,obj)
        
        assert idx is not None,"wrong index in add_object of Gene"

        self.intervals.update(obj.get_exons())
        self.intervals.merge_overlaps()

        self.start = min(self.start, obj.get_start())
        self.end = max(self.end, obj.get_end())

        return True
    
    def get_start(self) -> int:
        """Get the start coordinate of the Bundle.

        Returns:
            int: The start coordinate.
        """
        return self.start

    def get_end(self) -> int:
        """Get the end coordinate of the Bundle.

        Returns:
            int: The end coordinate.
        """
        return self.end

class Gene (Bundle):
    """
    Represents a specialization of Bundle where all objects have the same gene ID.
    """

    def __init__(self):
        super().__init__()
        self.gid = None

    def add_object(self,obj: Object) -> bool:
        """Add an object to the Gene.
        In addition to constraints from regular Bundle Gene asserts:
        1. all objects have the same gene_id

        If constrains fail - does not add object.

        Args:
            obj (Object): The object to add.

        Returns:
            bool: True if the object was successfully added, False otherwise.
        """
        obj_gid = obj.get_attr("gene_id")
        if obj_gid is None:
            return False
        
        if self.gid is None:
            self.gid = obj_gid

        if self.gid != obj_gid:
            return False
        
        status = super().add_object(obj)
        return status
    
    def get_gid(self) -> str:
        """Get the gene ID of the Gene.

        Returns:
            str: The gene ID.
        """
        return self.objects[0].get_gid()
    
class OverlapBundle(Bundle):
    """
    Represents a specialization of Bundle where every object added overlaps the bundle coordinates.
    It is up to the user to guarantee that transitively overlapping objects are added in the correct order
    """

    def __init__(self):
        super().__init__()

    def add_object(self,obj: Object) -> bool:
        """Add an object to the OverlapBundle.
        In addition to constraints from regular Bundle Gene asserts:
        1. every added object overlaps the current bundle coordinates.

        If constrains fail - does not add object.

        Args:
            obj (Object): The object to add.

        Returns:
            bool: True if the object was successfully added, False otherwise.
        """
        if self.seqid is None or self.overlaps(obj):
            status = super().add_object(obj)
            return status
        return False    

# largest specialization of the TXGroup which can yield all other types
class Transcriptome (TXGroup):
    """
    Represents the largest specialization of TXGroup that can hold and yield all other types of objects.
    Can be used to hold all objects together.
    Useful when loading unsorted files and need to make certain all child-parent relationships are recovered.
    """

    def __init__(self):
        super().__init__()

        self.genome = None
        
        self.donors = {}
        self.acceptors = {}

    def build_from_file(self,fname: str) -> None:
        """
        Build the Transcriptome from a file.
        Will try adding any object to the transcriptome groupping them by transcript id
        Genes and other collections of objects will be broken down and added as their constituent parts

        Args:
            fname (str): The file name to read the objects from.
        """
        treader = TReader(fname)
        for obj in treader.next_obj():
            self.add_object(obj)

        for obj in self.objects:
            if obj.get_type() == Types.Transcript:
                obj.finalize()

    def load_expression(self,fname: str) -> None:
        """
        Load expression data into the Transcriptome.
        Expression data is expected to be in the TSV format with column as follows:
        1. transcript_id
        2..N expression for each sample
        Each row is expected to contain the suame number of values

        Args:
            fname (str): The file name of the expression data in TSV format as described above.
        """
        assert os.path.exists(fname),"expression data file not found: "+fname

        num_samples = None

        with open(fname,"r") as inFP:
            for line in inFP:
                lcs = line.strip().split("\t")
                assert len(lcs)>1,"the expressiontable file format is expected to have more than 1 column: 1st columns is transcript ID and each successive column is expression in the sample"
                if num_samples is None:
                    num_samples = len(lcs)-1
                assert num_samples == len(lcs)-1,"number of samples is inconsistend in the expression file. Each row is expected to have the same number of values"
                
                # check if the line is header line
                has_tx_id_col = lcs[0]=="tx_id"
                numeric_values = [x.replace(".","",1).isnumeric() for x in lcs[1:]]

                if has_tx_id_col or not all(numeric_values):
                    continue

                cur_tid = lcs[0]
                idx = self.tid_map.get(cur_tid,None)

                if idx is not None:
                    for expstr in lcs[1:]:
                        exp = 0
                        try:
                            exp = float(expstr)
                        except:
                            pass
                        self.objects[idx].add_expression(exp)
                    
        
    def coordinate_sort(self):
        """
        Sort the objects in the Transcriptome based on coordinates.
        """
        self.objects.sort(key=lambda obj: obj.get_exons())
        self.reindex()

    def gid_sort(self):
        """
        Sort the objects in the Transcriptome based on gene ID.
        """
        self.objects.sort(key=lambda obj: obj.get_gid())
        self.reindex()

    def transcript_it(self):
        """
        Iterate over the Transcript objects in the Transcriptome.
        """
        for obj in self.object_it():
            if obj.obj_type == Types.Transcript:
                yield obj

    def gene_it(self) -> Gene:
        """
        Iterate over the Gene objects in the Transcriptome.
        Does not perform sorting of any kind - users are expected to guarantee themselves that transcriptome is appropriately sorted.
        """
        gene = Gene()
        for obj in self.object_it():
            if gene.is_empty() or gene.get_gid() == obj.get_gid():
                gene.add_object(obj)
            else:
                res_gene = copy.deepcopy(gene)
                gene = Gene()
                gene.add_object(obj)
                yield res_gene

        if not gene.is_empty():
            yield gene

    def bundle_it(self) -> OverlapBundle:
        """
        Creates and yields OverlapBundles.
        Does not perform sorting of any kind - users are expected to guarantee themselves that transcriptome is appropriately sorted.
        """
        bundle = OverlapBundle()
        for obj in self.object_it():
            if bundle.is_empty() or bundle.overlaps(obj):
                bundle.add_object(obj)
            else:
                res_bundle = copy.deepcopy(bundle)
                bundle = OverlapBundle()
                bundle.add_object(obj)
                yield res_bundle
                

        if not bundle.is_empty():
            yield bundle

    def load_genome(self, fasta_fname):
        self.genome = Fasta(fasta_fname)
        
    def extract_introns(self):
        for tx in self.transcript_it():
            for it in tx.introns_it():
                self.donors.setdefault(it[0],{"count":0,"seq":None})
                self.donors[it[0]]["count"]+=1
                
                self.acceptors.setdefault(it[0],{"count":0,"seq":None})
                self.acceptors[it[0]]["count"]+=1

    def intervals(self):
        """
        Returns:
            IntervalTree: A tree of all exons in the Transcriptome.
        """
        it = IntervalTree()
        for obj in self.transcript_it():
            it.update(obj.get_exons())
        it.merge_overlaps()
        return it