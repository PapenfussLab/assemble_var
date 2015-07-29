"""
fasta module
"""

from mungoCore import *
from useful import smartopen, progressMessage


class Motif(AbstractFeature):
    attributes = ['factor','start','end','strand','seq','score']
    converters = [
        ('factor', None),
        ('start', int),
        ('end', int),
        ('strand', None),
        ('seq', None),
        ('score', float)
    ]
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        super(Motif,self).__init__(*args,**kw)


class PossumOutputReader(AbstractDataReader):
    """Simple class for reading possum output files"""
    
    def __init__(self, iFileHandle, **kw):
        """
        @param iFileHandle: Fasta file name or object
        """
        AbstractDataReader.__init__(self, iFileHandle)
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def _generator(self):
        header = ''
        motifs = []
        isStart = True
        for line in self.iFile:
            line = line.strip()
            if line:
                if isStart and line[0]!='>': 
                    continue
                else:
                    isStart = False
                
                if line[0]=='>':
                    if motifs:
                        yield header, motifs
                    header = line[1:]
                    motifs = []
                else:
                    tokens = line.strip().split()
                    del tokens[2]
                    motif = Motif(tokens)
                    motifs.append(motif)
        if header and motifs:
            yield header, ''.join(motifs)
