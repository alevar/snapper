from typing import Iterator

from ..classes.transcript import GTFObjectFactory, Object
from ..utils.common import *

class TReader:
    """
    A class to read and validate general structure of GTF/GFF files

    ...

    Attributes
    ----------
    fname : str
        filename of the GTF/GFF file to be read
    fp : File
        file handler opened automatically by the TReader when opening the input files
    gff : bool
        set to True if the file is detected to be in the GFF format. Set to false if in GTF
    comments : list
        list of all comments found in the file is the same order as they were read in the input file.

    Methods
    -------
    _set_file(fname="") : None
        Internal function that sets fname for a TReader object and runs format detection setting the gff attribute accordingly
    _is_gff() : bool
        runs automatic detection of the input file format
    is_gff() : bool
        returns true if the file is detected to be in GFF format
    next_obj() : Iterator[Object]
        advances the file reader by single line skipping comments. Uses a factory to construct appropriate objects
    """
    def __init__(self,fname:str=None):
        self.fname = ""
        self.fp = None
        self.gff = False
        self.comments = []
        
        if not fname is None:
            self._set_file(fname)

    def __del__(self):
        if not self.fp is None:
            self.fp.close()

    def _set_file(self,fname:str):
        self.fname = fname
        if self.fp is not None:
            self.fp.close()

        self.gff = self._is_gff()
        if self.gff is None:
            raise Exception("File is not a valid GTF or GFF file: "+fname)
        self.fp = open(self.fname,"r")

    # function opens the file, peaks inside, and returns True if the file is in GTF format, False if it is in GFF format and None if it is not a valid file
    # will always close the file and re-open it
    def _is_gff(self) -> bool:
        gff = None
        if not self.fp is None:
            self.fp.close()

        self.fp = open(self.fname,"r")

        for line in self.fp:
            if line.startswith("#"):
                continue
            lcs = line.strip().split("\t")
            if len(lcs) > 9:
                gff = None
                break
            
            if lcs[2] not in ["transcript","exon","CDS"]:
                gff = None
                break

            if lcs[2] == "transcript":
                if lcs[8].startswith("ID="):
                    gff = True
                    break
                elif lcs[8].startswith("transcript_id"):
                    gff = False
                    break
                else:
                    gff = None
                    break
            else:
                continue
        
        self.fp.close()
        self.fp = open(self.fname,"r")
        return self.is_gff()
    
    def is_gff(self):
        return self.gff

    def next_obj(self) -> Iterator[Object]:
        for line in self.fp:
            if line.startswith('#'):
                self.comments.append(line)
                continue

            obj = GTFObjectFactory.create(line)
            if not obj is None:
                yield obj