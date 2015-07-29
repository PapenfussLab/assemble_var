"""
fastq module
"""

from mungoCore import *
from useful import smartopen


def FastqFile(fileHandle, mode='r', **kw):
    """Factory function for Reader and Writer classes.
    
    @param iFileHandle: Fasta file name or object
    @keyword mode: read(r), write(w) or append(a)
    @keyword indexed: Index the fasta file; reading only (default False)
    """
    if mode=='r':
        return FastqReader(fileHandle, **kw)
    elif mode=="w":
        return FastqWriter(fileHandle, **kw)
    else:
        raise NotImplemented


class FastqReader(AbstractDataReader):
    """Simple class for reading fastq files"""
    
    def __init__(self, iFileHandle, **kw):
        """
        @param iFileHandle: Fasta file name or object
        """
        AbstractDataReader.__init__(self, iFileHandle)
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def __iter__(self):
        self._iter = self._generator()
        return self
    
    def next(self):
        for h,s,q in self._iter:
            return h,s,q
        raise StopIteration
    
    def _generator(self):
        header = ""
        qheader = ""
        seq = []
        qual = []
        state = "B"
        for line in self.iFile:
            line = line.strip()
            if line:
                if state=="Q" and line[0]=="@":  # Internal read header
                    yield header, "".join(seq), "".join(qual)
                    state = "S"
                    header = line[1:]
                    seq = []
                    qual = []
                elif state=="S" and line[0]=="+":  # Internal quality header
                    state = "Q"
                    header = line[1:]
                elif state=="B" and line[0]=="@": # Starting read header
                    state = "S"
                    header = line[1:]
                    seq = []
                    qual = []
                else: # Internal read or quality sequence
                    if state=="S":
                        seq.append(line)
                    elif state=="Q":
                        qual.append(line)
        
        if state=="Q":
            yield header, "".join(seq), "".join(qual)


class FastqWriter(AbstractDataFile):
    """Simple class for writing fastq files"""

    def __init__(self, fileHandle, mode='w', width=60, blockSize=None, **kw):
        """
        @param iFileHandle: Output file or name
        @keyword mode: File mode - write(w) or append(a)
        """
        assert mode in ('w', 'a')
        self.iFile = smartopen(fileHandle, mode)
        self.iFilename = self.iFile.name
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def write(self, header, seq, qseq):
        """Write a fasta entry.
        
        @param header: Fasta header
        @param seq: Sequence
        @param qseq: Quality sequence
        """
        print >> self.iFile, "@%s" % header
        print >> self.iFile, seq
        print >> self.iFile, "+%s" % header
        print >> self.iFile, qseq
    
    def __call__(self, header, seq, qseq):
        """Call interface to write."""
        self.write(header, seq, qseq)
