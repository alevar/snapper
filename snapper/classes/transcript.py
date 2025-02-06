from intervaltree import Interval, IntervalTree
from ..utils.common import *
from typing import Tuple,List, Callable
import copy

class Object:
    def __init__(self):
        self.start = float("inf") # inclusive
        self.end = 0 # non-inclusive

        self.source = "ORFclust"
        self.obj_type = None
        
        self.seqid = None
        self.strand = None

        self.attrs = dict()
        self.expression = list()
        
        self.data = None # whatever data user wants to store in the object

    def clear(self):
        """
        Clear the attributes of the object.
        """
        self.seqid = None
        self.strand = None
        self.source = "ORFclust"
        self.obj_type = None
        self.start = float("inf")
        self.end = 0
        self.attrs.clear()
        self.expression = list()
        self.data = None

    def to_transcript(self) -> 'Transcript':
        """
        Convert the object to a Transcript object.

        Returns:
            Transcript: A Transcript object.
        """
        tx = Transcript()
        tx.set_seqid(self.seqid)
        tx.set_strand(self.strand)
        tx.set_source(self.source)
        tx.set_type(Types.Transcript)
        tx.set_start(self.start)
        tx.set_end(self.end)
        tx.set_attributes(self.attrs)
        tx.set_tid(self.attrs["transcript_id"])
        tx.set_gid(self.attrs.get("gene_id", None))
        tx.set_expression(self.expression)
        tx.set_data(self.data)
        return tx
    
    def to_exon(self) -> 'Exon':
        """
        Convert the object to an Exon object.

        Returns:
            Exon: An Exon object.
        """
        exon = Exon()
        exon.set_seqid(self.seqid)
        exon.set_strand(self.strand)
        exon.set_source(self.source)
        exon.set_type(Types.Exon)
        exon.set_start(self.start)
        exon.set_end(self.end)
        exon.set_attributes(self.attrs)
        exon.set_tid(self.attrs["transcript_id"])
        exon.set_gid(self.attrs.get("gene_id", None))
        exon.set_expression(self.expression)
        exon.set_data(self.data)
        return exon

    def is_empty(self) -> bool:
        """
        Check if the object is empty.

        Returns:
            bool: True if the object is empty, False otherwise.
        """
        return self.start == float("inf")

    def overlaps(self, obj) -> bool:
        """
        Check if the object overlaps with another object.

        Args:
            obj (Object): Another object to check for overlap.

        Returns:
            bool: True if the objects overlap, False otherwise.
        """
        return (self.get_start()<=obj.get_end() and self.get_end()>=obj.get_start())
    
    def set_seqid(self, seqid: str) -> None:
        """
        Set the sequence ID of the object.

        Args:
            seqid (str): The sequence ID.

        Returns:
            None
        """
        self.seqid = seqid

    def set_strand(self, strand: str) -> None:
        """
        Set the strand of the object.

        Args:
            strand (str): The strand.

        Returns:
            None
        """
        self.strand = strand

    def set_source(self, source: str) -> None:
        """
        Set the source of the object.

        Args:
            source (str): The source.

        Returns:
            None
        """
        self.source = source

    def set_type(self, obj_type: str) -> None:
        """
        Set the type of the object.

        Args:
            obj_type (str): The object type.

        Returns:
            None
        """
        self.obj_type = obj_type

    def set_start(self, start: int):
        """
        Set the start position of the object.

        Args:
            start (int): The start position.

        Returns:
            None
        """
        self.start = start

    def set_end(self, end: int):
        """
        Set the end position of the object.

        Args:
            end (int): The end position.

        Returns:
            None
        """
        self.end = end

    def set_attributes(self, attrs: dict) -> None:
        """
        Set the attributes of the object.

        Args:
            attrs (dict): A dictionary of attributes.

        Returns:
            None
        """
        self.attrs = attrs

    def add_attribute(self, k: str, v: str, replace: bool=False, append: bool=False) -> None: # replace - if key exists - replace it with new. append - if key exists - add new values to the list of values. If both enabled - will replace
        """
        Add or update an attribute to the object.

        Args:
            k (str): The attribute key.
            v (str): The attribute value.
            replace (bool, optional): If True and the key exists, replace the value. Defaults to False.
            append (bool, optional): If True and the key exists, append the value to the list. Defaults to False.

        Returns:
            None
        """
        if append and not replace:
            if k in self.attrs:
                old_val = self.attrs[k]
                if not old_val == v:
                    self.attrs[k] = old_val+", "+v # TODO: should we replace with a list instead?
        else:
            self.attrs.setdefault(k,v)
            if replace:
                self.attrs[k] = v

    def add_expression(self, exp: float) -> None:
        """
        Add an expression value to the object.

        Args:
            exp (float): The expression value.

        Returns:
            None
        """
        self.expression.append(exp)

    def set_expression(self, exps: list) -> None:
        """
        Set the expression values of the object.

        Args:
            exp (list): A list of expression values.

        Returns:
            None
        """
        self.expression = exps
        
    def set_data(self, data):
        """
        Set the data of the object.

        Args:
            data: The data.

        Returns:
            None
        """
        self.data = data

    def get_source(self) -> str:
        """
        Get the source of the object.

        Returns:
            str: The source.
        """
        return self.source
    
    def get_type(self) -> str:
        """
        Get the type of the object.

        Returns:
            str: The object type.
        """
        return self.obj_type

    def get_seqid(self) -> str:
        """
        Get the sequence ID of the object.

        Returns:
            str: The sequence ID.
        """
        return self.seqid

    def get_strand(self) -> str:
        """
        Get the strand of the object.

        Returns:
            str: The strand.
        """
        return self.strand

    def get_start(self) -> int:
        """
        Get the start position of the object.

        Returns:
            int: The start position.
        """
        return self.start

    def get_end(self) -> int:
        """
        Get the end position of the object.

        Returns:
            int: The end position.
        """
        return self.end

    def get_attr(self, attr: str) -> str:
        """
        Get the value of a specific attribute.

        Args:
            attr (str): The attribute name.

        Returns:
            str: The attribute value, or None if the attribute does not exist.
        """
        return self.attrs.get(attr, None)

    def get_attributes(self) -> dict:
        """
        Get all attributes of the object.

        Returns:
            dict: A dictionary of attributes.
        """
        return self.attrs

    def get_expression(self,func: Callable=None) -> list:
        """
        Get the expression values of the object.

        Args:
            func (Callable, optional): A function to apply to the expression values. Defaults to None.

        Returns:
            list: A list of expression values.
        """
        return func(self.expression) if func is not None else self.expression

    def add_line(self,gtf_line:str,store_all_attributes=False) -> bool:
        """
        Add information from a GTF line to the object.

        Args:
            gtf_line (str): The GTF line.
            store_all_attributes (bool, optional): If True, store all attributes including exon and CDS attributes.
                                                  Defaults to False.

        Returns:
            bool: True if the line was successfully added, False otherwise. Returns None in case of errors
        """
        lcs = gtf_line.strip().split("\t")
        if not len(lcs) == 9:
            return None

        self.attrs = extract_attributes(lcs[8])

        lstart = int(lcs[3])
        lend = int(lcs[4])+1 # +1 here because intervaltree inclusive of the lower limit, but non-inclusive of the upper limit
        
        self.start = min(self.start,lstart)
        self.end = max(self.end,lend)
        self.seqid = lcs[0]
        self.strand = lcs[6]

        if lcs[2] == "transcript" or lcs[2] == "mRNA":
            self.obj_type = Types.Transcript
            self.attrs = rename_attributes(self.attrs,{"ID":"transcript_id","Parent":"gene_id"})
        elif lcs[2] == "exon":
            self.obj_type = Types.Exon
            self.attrs = rename_attributes(self.attrs,{"Parent":"transcript_id"})
        elif lcs[2] == "CDS":
            self.obj_type = Types.CDS
            self.attrs = rename_attributes(self.attrs,{"Parent":"transcript_id"})
        else:
            self.obj_type = Types.Other
        
        return True
    
    def _getattr(self, name):
        """
        Allows queries of both the object and the attributes parsed from the file
        """
        try:
            return getattr(self,name)
        except:
            if name in self.attrs:
                return self.attrs[name]
            else:
                raise AttributeError(f"'Object' object has no attribute '{name}'")
            
    def _getattrs(self,attrs) -> Tuple:
        """
        Able to return one or more attributes at once. Used when sorting or groupping by multiple keys.
        """
        return tuple([self._getattr(attr) for attr in attrs] if isinstance(attrs, list) else [self._getattr(attrs)])
    
    def copy(self):
        """
        Create a copy of the object.

        Returns:
            Object: A new instance of the Object class with the same attribute values.
        """
        obj = Object()
        obj.set_seqid(self.seqid)
        obj.set_source(self.source)
        obj.set_strand(self.strand)
        obj.set_type(self.obj_type)
        obj.set_attributes(self.attrs)
        obj.set_start(self.start)
        obj.set_end(self.end)
        obj.expression = [x for x in self.expression]
        return obj
    
    def len(self):
        """
        Get the length of the object.

        Returns:
            int: The length of the object.
        """
        return self.end-self.start
    
    def get_sequence(self,genome):
        sequence = genome[self.seqid][self.start-1:self.end-1].seq
        if self.strand == '-':
            sequence = reverse_complement(sequence)
        return sequence

    def genome_coordinate(self, pos: int) -> int:
        """
        Convert a relative position to a genomic coordinate.

        Args:
            pos (int): The relative position.

        Returns:
            int: The genomic coordinate.
        """
        return self.start+pos
    
    def merge(self,obj: 'Object'):
        """
        Merge another object into the object.
        keeps IDs and data of the self object, only modifies coordinates and attributes

        Args:
            obj (Object): The object to merge.

        Returns:
            None
        """
        self.start = min(self.start,obj.get_start())
        self.end = max(self.end,obj.get_end())

class Transcript (Object):
    """
    A class representing a transcript.

    Attributes:
        obj_type (Types): The type of the object (Transcript).
        exons (IntervalTree): An interval tree containing exons.
        cds (IntervalTree): An interval tree containing CDS regions.

    Inherits from:
        Object

    """
    def __init__(self, obj: Object=None):
        super().__init__()
        self.tid = None
        self.gid = None
        self.obj_type = Types.Transcript
        self.exons = IntervalTree()
        self.cds = IntervalTree()

        if obj is not None:
            self.__dict__ = obj.to_transcript().__dict__.copy()

    def merge(self,obj: Object, merge_cds:bool=True):
        """
        Merge object into the transcript. Ensured transcript_ids match. Can merge Exons, CDS, and Transcripts into the Transcript.
        
        Args:
            obj (Object): The Object to merge.
            merge_cds (bool, optional): Whether to merge CDS regions. Defaults to True.
        """
        # assert self.tid is None or self.tid == obj.get_attr("transcript_id"), "Transcript IDs do not match"
        assert self.seqid is None or self.seqid == obj.get_seqid(), "Sequence IDs do not match"
        assert self.strand is None or self.strand == obj.get_strand(), "Strands do not match"
        
        # # if gene_id is not set, but the transcript_id is set and matches a known transcript_id - match gene_id
        # if self.gid is None and obj.get_attr("gene_id") is not None:
        #     self.gid = obj.get_attr("gene_id")
        # elif self.gid is not None and obj.get_attr("gene_id") is None:
        #     obj.add_attribute("gene_id",self.gid)
        # else:
        #     assert self.gid == obj.get_attr("gene_id"), "Gene IDs do not match"
            
        for k,v in obj.get_attributes().items():
            if k in self.attrs:
                continue
            self.add_attribute(k,v)
        if obj.get_type() == Types.Exon:
            self.add_exon(obj)
        elif obj.get_type() == Types.CDS:
            self.add_cds(obj)
        elif obj.get_type() == Types.Transcript:
            for exon in obj.get_exons():
                self.add_exon(exon[2])
            self.merge_exons()
            if merge_cds:
                for cds in obj.get_cds():
                    self.add_cds(cds[2])
                self.merge_cds()

    def own(self, obj: Object):
        """
        Converts attributes from the object to the current object
        
        Args:
            obj (Object): The Object to own.
        """
        obj.set_tid(self.get_tid())
        obj.set_gid(self.get_gid())

    def add_exon(self,obj: Object) -> bool: # returns True if sucessfully added
        """
        Add an Exon to the Transcript. Objects are sorted accordingly to update coordinates, boundaries, exon and cds chains

        Args:
            obj (Object): The Object to add.
            
        Returns:
            bool: True if the Object was successfully added, False otherwise.

        """
        if self.strand != obj.get_strand() or self.seqid != obj.get_seqid():
            return False
        # if self.tid != obj.get_attr("transcript_id"):
        #     return False
        
        exon = obj.to_exon()
        self.own(exon)
        self.start = min(self.start,exon.get_start())
        self.end = max(self.end,exon.get_end())
        self.exons.addi(exon.get_start(),exon.get_end(),exon)

        return True
    
    def add_cds(self,obj: Object) -> bool: # returns True if sucessfully added
        """
        Add an CDS to the Transcript. Objects are sorted accordingly to update coordinates, boundaries, exon and cds chains

        Args:
            obj (Object): The Object to add.

        Returns:
            bool: True if the Object was successfully added, False otherwise.

        """
        if self.strand != obj.get_strand() or self.seqid != obj.get_seqid():
            return False
        if self.tid != obj.get_attr("transcript_id"):
            return False
        
        cds = CDS(obj.to_exon())
        assert cds.get_start() >= self.start and cds.get_end() <= self.end,"CDS out of transcript boundaries"
        self.cds.addi(cds.get_start(),cds.get_end(),cds)

        return True

    def finalize(self, extend: bool=False) -> None:
        """
        Make sure that the transcript is complete.

        If extend is enabled, the 3' and 5' exons will be extended to match the start and end of the transcript.
        # consider that only a transcript line has been added - this function will create required exons and intervals
        # if exon boundaries and transcript boundaries are out of sync (transcript line describes longer interval than exon chain) - they will be adjusted accordingly

        Args:
            extend (bool, optional): Whether to extend the exons. Defaults to False.

        Returns:
            None

        """
        assert self.tid == self.attrs["transcript_id"],"transcript_id missing - can not assign"
        if len(self.exons)>0:
            exon_start = sorted(self.exons)[0][0]
            exon_end = sorted(self.exons)[-1][1]
            if extend:
                self.start = min(self.start,exon_start)
                self.end = max(self.end,exon_end)
                self.exons.smallest = self.start # TODO: implement
                self.exons.largest = self.end # TODO: implement
            else:
                self.start = exon_start
                self.end = exon_end
        else:
            obj = Object()
            obj.set_seqid(self.seqid)
            obj.set_strand(self.strand)
            obj.set_start(self.start)
            obj.set_end(self.end)
            obj.set_type(Types.Exon)
            obj.set_attributes({"transcript_id":self.tid})
            exon = Exon(obj)
            self.exons.addi(self.start,self.end,exon)

        if len(self.cds)>0:
            # make sure the CDS (if added) fits within exons
            assert len(self.exons[self.cds.begin()])>0 and len(self.exons[self.cds.end()-1])>0,"invalid CDS detected in transcript when finalizing: "+self.tid

            self.assign_phase()

    @staticmethod
    def _intervaltree_data_reducer(how: str = "union") -> Callable:
        """Outer function that takes additional arguments."""
        def data_reducer(obj1, obj2):
            """Inner function that does the merging."""
            if how == "union":
                obj1.merge(obj2)
                return obj1
            elif how == "longest":
                return obj1 if obj1.len() > obj2.len() else obj2
            elif how == "first":
                # return the object with the smallest start position and if start is the same, the one with the smallest end position
                return obj1 if obj1.get_start() < obj2.get_start() or (obj1.get_start() == obj2.get_start() and obj1.get_end() < obj2.get_end()) else obj2
            else:
                raise Exception("Invalid data reducer")
        
        return data_reducer

    # @staticmethod
    # def _intervaltree_data_reducer(obj1, obj2, how="union"):
    #     if how == "union":
    #         obj1.merge(obj2)
    #         return obj1
    #     elif how == "longest":
    #         return obj1 if obj1.len() > obj2.len() else obj2
    #     elif how == "first":
    #         return obj1 if obj1.get_start() < obj2.get_start() else obj2
    #     else:
    #         raise Exception("Invalid data reducer")

    def merge_exons(self,how="union") -> None:
        """
        Merge overlapping exons.

        Returns:
            None

        """
        reducer = Transcript._intervaltree_data_reducer(how)
        self.exons.merge_overlaps(data_reducer = reducer)
        
    def merge_cds(self,how="union") -> None:
        """
        Merge overlapping CDS regions.

        Returns:
            None

        """
        reducer = Transcript._intervaltree_data_reducer(how)
        self.cds.merge_overlaps(data_reducer = reducer)

    def assign_phase(self, start_phase: int=0) -> None:
        """
        Assign the phase to the CDS elements.

        Args:
            start_phase (int, optional): The phase to start with. Defaults to 0.
        """

        if self.strand == "-":
            cdsacc = start_phase
            for c in sorted(self.cds)[::-1]:
                c[2].set_phase((3-cdsacc%3)%3)
                cdsacc += (c[2].get_end()-1) - c[2].get_start()+1
        else:
            cdsacc = start_phase
            for c in sorted(self.cds):
                c[2].set_phase((3-cdsacc%3)%3)
                cdsacc += (c[2].get_end()-1) - c[2].get_start()+1

    def to_exon(self) -> 'Exon':
        """
        Convert the Transcript to an Exon.

        Returns:
            Exon: The Exon.

        """
        obj = Object()
        obj.set_seqid(self.seqid)
        obj.set_strand(self.strand)
        obj.set_start(self.start)
        obj.set_end(self.end)
        obj.set_type(Types.Exon)
        obj.set_attributes({"transcript_id":self.tid})
        exon = Exon(obj)
        return exon
    
    def to_transcript(self) -> 'Transcript':
        """
        Convert the Transcript to a Transcript.

        Returns:
            Transcript: The Transcript.

        """
        return self
    
    def clear(self):
        """
        Clear the Transcript object.

        Returns:
            None

        """
        super().clear()

        self.exons = IntervalTree()
        self.cds = IntervalTree()
        self.tid = None
        self.gid = None
        self.obj_type = Types.Transcript

    def set_exons(self, exons: list[tuple[int, int]]) -> None:
        """
        Set the exons of the Transcript.
        Objects will be copied with their attributes including transcript and gene ID.
        Up to the user to ensure attributes are handled correctly.

        Args:
            exons (list[tuple[int, int]]): A list of exon intervals represented as tuples.

        Returns:
            None

        """
        for e in exons:
            exon_obj = Exon()
            if len(e)==3 and type(e[2])==Exon: # has Exon object in the tuple
                exon_obj = copy.deepcopy(e[2])
            elif len(e)>=2:
                exon_obj.set_tid(self.tid)
                exon_obj.set_seqid(self.seqid)
                exon_obj.set_strand(self.strand)
                exon_obj.set_start(e[0])
                exon_obj.set_end(e[1])
                exon_obj.set_attributes(e[2].get_attributes())
            else:
                raise Exception("invalid Exon tuple: "+str(e))
            
            self.exons.addi(e[0],e[1],exon_obj)

    def set_cds(self, cds: list[tuple[int, int]]) -> None:
        """
        Set the CDS regions of the Transcript.
        Objects will be copied with their attributes including transcript and gene ID.
        Up to the user to ensure attributes are handled correctly.

        Args:
            cds (list[tuple[int, int]]): A list of CDS intervals represented as tuples.

        Returns:
            None

        """
        for c in cds:
            cds_obj = CDS()
            if len(c)==3 and type(c[2])==CDS: # has CDS object in the tuple
                cds_obj = copy.deepcopy(c[2])
            elif len(c)>=2:
                cds_obj.set_tid(self.tid)
                cds_obj.set_seqid(self.seqid)
                cds_obj.set_strand(self.strand)
                cds_obj.set_start(c[0])
                cds_obj.set_end(c[1])
                cds_obj.set_phase(0)
                cds_obj.set_attributes(c[2].get_attributes())
            else:
                raise Exception("invalid CDS tuple: "+str(c))
            
            self.cds.addi(c[0],c[1],cds_obj)

            # TODO: reassign phase when finalizing transcript

    def clear_exons(self) -> None:
        """
        Clear the exons of the Transcript.

        Returns:
            None

        """
        self.exons = IntervalTree()
    
    def clear_cds(self) -> None:
        """
        Clear the CDS regions of the Transcript.

        Returns:
            None

        """
        self.cds = IntervalTree()

    def set_tid(self, tid: str) -> None:
        """
        Set the transcript ID of the Transcript.

        Args:
            tid (str): The transcript ID.

        Returns:
            None

        """
        self.tid = tid
        self.add_attribute("transcript_id",tid,replace=True)

    def set_gid(self, gid: str) -> None:
        """
        Set the gene ID of the Transcript.

        Args:
            gid (str): The gene ID.

        Returns:
            None

        """
        self.gid = gid
        self.add_attribute("gene_id",gid,replace=True)

    def nume(self):
        """
        Get the number of exons in the Transcript.

        Returns:
            int: The number of exons.

        """
        return len(self.exons)

    def numc(self):
        """
        Get the number of CDS regions in the Transcript.

        Returns:
            int: The number of CDS regions.

        """
        return len(self.cds)

    def has_cds(self):
        """
        Check if the Transcript has CDS regions.

        Returns:
            bool: True if the Transcript has CDS regions, False otherwise.

        """
        return len(self.cds) > 0

    def introns_it(self):
        """Yield all introns of the transcript.

        Intervals are also inclusive of the lower bound but non-inclusive of the upper bound.
        [first position of the intron, last position of the intron)

        Yields:
            Interval: An Interval object representing an intron.

        """
        if len(self.exons) > 1:
            prev_exon = None
            for e in sorted(self.exons):
                if prev_exon is not None:
                    yield Interval(prev_exon[1], e[0])
                prev_exon = e

    def get_exons(self):
        """Get the sorted list of exons in the Transcript.

        Returns:
            list: A list of exon intervals represented as tuples.

        """
        return sorted(self.exons)

    def get_cds(self):
        """Get the sorted list of CDS regions in the Transcript.

        Returns:
            list: A list of CDS intervals represented as tuples.

        """
        return sorted(self.cds)

    def get_tid(self):
        """Get the transcript ID of the Transcript.

        Returns:
            str: The transcript ID.

        """
        return self.tid

    def get_gid(self):
        """Get the gene ID of the Transcript.

        Returns:
            str: The gene ID.

        """
        return self.gid
    
    def get_cstart(self):
        """Get the start position of the first CDS region in the Transcript.

        Returns:
            int: The start position of the first CDS region.

        """
        return sorted(self.cds)[0][0]

    def get_cend(self):
        """Get the end position of the last CDS region in the Transcript.

        Returns:
            int: The end position of the last CDS region.

        """
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        """Convert the Transcript to GTF format.

        Returns:
            str: The Transcript object represented in GTF format.

        """
        res =  self.seqid+"\t"+\
               self.source+"\t"+\
               Types.type2str(self.obj_type) +"\t"+\
               str(self.start)+"\t"+\
               str(self.end-1)+"\t"+\
               "."+"\t"+\
               self.strand+"\t"+\
               "."+"\t"+\
               to_attribute_string(self.attrs,False,"transcript")+"\n"
        
        for e in sorted(self.exons):
            res += e[2].to_gtf() + "\n"
            
        for c in sorted(self.cds):
            res += c[2].to_gtf() + "\n"

        return res.rstrip("\n")

    def to_gff(self):
        """Convert the Transcript to GFF format.

        Returns:
            str: The Transcript object represented in GFF format.

        """
        res =  self.seqid+"\t"+\
               self.source+"\t"+\
               Types.type2str(self.obj_type) +"\t"+\
               str(self.start)+"\t"+\
               str(self.end-1)+"\t"+\
               "."+"\t"+\
               self.strand+"\t"+\
               "."+"\t"+\
               to_attribute_string(self.attrs,True,"transcript")+"\n"
        
        for e in sorted(self.exons):
            res += e[2].to_gff() + "\n"
            
        for c in sorted(self.cds):
            res += c[2].to_gff() + "\n"

        return res.rstrip("\n")

    __repr__ = to_gtf

    __str__ = __repr__

    def copy(self):
        tx = Transcript()
        tx.set_tid(self.tid)
        tx.set_gid(self.gid)
        tx.set_exons(self.exons)
        tx.set_cds(self.cds)
        return tx
    
    def elen(self) -> int:
        """Get effective length of the transcript in base pairs.

        Returns:
            int: The length of the transcript.

        """
        return sum([e[2].len() for e in self.exons])
    
    def clen(self) -> int:
        """Get coding length of the transcript in base pairs.

        Returns:
            int: The length of the ORF.

        """
        return sum([c[2].len() for c in self.cds])
    
    def get_sequence(self, genome, use_cds: bool=False) -> str:
        sequence = ""
        # Check if the strand is minus
        objects = sorted(self.cds if use_cds else self.exons, reverse=self.strand == '-')
        for o in objects:
            sequence += o[2].get_sequence(genome)
        return sequence

    def genome_coordinate(self, pos: int) -> int:
        """Convert a relative position to a genomic coordinate.

        Args:
            pos (int): The relative position.

        Returns:
            int: The genomic coordinate.

        """
        cur_len = 0
        for e in sorted(self.exons):
            if cur_len + e[2].len() > pos:
                return e[2].genome_coordinate(pos-cur_len)
            cur_len += e[2].len()
        return None

    def get_chain(self, use_cds: bool=False) -> List[int]:
        """Get the chain of the transcript.

        Args:
            use_cds (bool, optional): If True, use the CDS chain. Defaults to False.

        Returns:
            List[int]: A list of the exon coordinates.

        """
        chain = []
        objects = sorted(self.cds if use_cds else self.exons, reverse=self.strand == '-')
        for o in objects:
            chain.append([o[2].get_start(), o[2].get_end()])
        return chain
    
class Exon(Object):
    """
    A class representing an exon.

    Attributes:
        obj_type (Types): The type of the object (Exon).

    Inherits from:
        Object

    """

    def __init__(self,obj: Object=None):
        """Initialize an Exon object.

        Args:
            obj (Object, optional): An Object from which to initialize the Exon. Defaults to None.

        Raises:
            Exception: If a wrong object type is passed to the constructor.

        """
        super().__init__()
        self.obj_type = Types.Exon
        self.tid = None
        self.gid = None
        if obj is not None:
            self.__dict__ = obj.to_exon().__dict__.copy()

    def set_tid(self, tid: str) -> None:
        """Set the transcript ID.

        Args:
            tid (str): The transcript ID.

        """
        self.tid = tid
        self.add_attribute("transcript_id",tid,replace=True)

    def set_gid(self, gid: str) -> None:
        """Set the gene ID.

        Args:
            gid (str): The gene ID.

        """
        self.gid = gid
        self.add_attribute("gene_id",gid,replace=True)

    def get_tid(self) -> str:
        """Get the transcript ID.

        Returns:
            str: The transcript ID.

        """
        return self.tid

    def get_gid(self) -> str:
        """Get the gene ID.

        Returns:
            str: The gene ID.

        """
        return self.gid

    def get_cstart(self):
        """Get the start position of the first CDS region in the Exon.

        Returns:
            int: The start position of the first CDS region.

        """
        return sorted(self.cds)[0][0]

    def get_cend(self):
        """Get the end position of the last CDS region in the Exon.

        Returns:
            int: The end position of the last CDS region.

        """
        return sorted(self.cds)[-1][1]

    def to_gtf(self):
        """Convert the Exon to GTF format.

        Returns:
            str: The Exon object represented in GTF format without new line at the end.

        """
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                to_attribute_string(self.attrs,False,"exon")
        return res
    
    def to_gff(self):
        """Convert the Exon to GFF format.

        Returns:
            str: The Exon object represented in GFF format  without new line at the end.

        """
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                "."+"\t"+\
                to_attribute_string(self.attrs,True,"exon")
        return res
    
    __repr__ = to_gtf

    __str__ = __repr__

    def __eq__(self, other: Object) -> bool:
        '''
        Check if two objects are equal. Only consider coordinates and not attributes or tid
        '''
        return isinstance(other, Object) and self.seqid==other.seqid and self.start==other.start and self.end==other.end and self.strand==other.strand
    

class CDS(Exon):

    def __init__(self,obj: Object=None):
        """Initialize a CDS object.

        Args:
            obj (Object, optional): An Object from which to initialize the CDS. Defaults to None.

        Raises:
            Exception: If a wrong object type is passed to the constructor.

        """
        super().__init__(obj)
        self.obj_type = Types.CDS
        self.phase = 0
    
    def set_phase(self,phase):
        self.phase = phase

    def get_phase(self):
        return self.phase
    
    def to_gtf(self):
        """Convert the CDS to GTF format.

        Returns:
            str: The CDS object represented in GTF format without new line at the end.

        """
        res  =  self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                str(self.phase)+"\t"+\
                to_attribute_string(self.attrs,False,"exon")
        return res
    
    def to_gff(self):
        """Convert the CDS to GFF format.

        Returns:
            str: The CDS object represented in GFF format  without new line at the end.

        """
        res = self.seqid+"\t"+\
                self.source+"\t"+\
                Types.type2str(self.obj_type) +"\t"+\
                str(self.start)+"\t"+\
                str(self.end-1)+"\t"+\
                "."+"\t"+\
                self.strand+"\t"+\
                str(self.phase)+"\t"+\
                to_attribute_string(self.attrs,True,"exon")
        return res
    
    __repr__ = to_gtf

    __str__ = __repr__

class GTFObjectFactory:
    """
    Object Factory Class.
    Automatically generates objects of correct type from GTF/GFF lines

    """
    @staticmethod
    def create(line):
        """Create an object based on the GTF line.

        Args:
            line (str): A line from the GTF file.

        Returns:
            Object: The created object based on the GTF line.

        """
        obj = Object()
        obj.add_line(line)
        if obj.get_type() == Types.Transcript or obj.get_type() == Types.MRNA:
            return obj.to_transcript()
        elif obj.get_type() == Types.Exon:
            return obj.to_exon()
        elif obj.get_type() == Types.CDS:
            return CDS(obj)
        else:
            return obj
        
    def __eq__(self, other: Object) -> bool:
        '''
        Check if two objects are equal. Only consider coordinates and not attributes or tid
        '''
        return super().__eq__() and self.phase==other.phase
    


# every GTF object has some shared properties
# every object has coordinate range seqid, strand, start and end.
# every object can not span multiple chromosomes and seqids
# every object has source (otherwise set as default)
# every object can (but does not have to) have attributes