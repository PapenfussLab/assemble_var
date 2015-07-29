"""
blat module
"""

from mungoCore import *
import sequence
from useful import smartopen, DSU


class Chain(AbstractFeature):
    """Blat chain feature"""
    
    attributes = [
        'match','mismatch','repMatch','Ns',
        'qGapCount','qGapBases',
        'tGapCount','tGapBases','strand',
        'qName','qSize','qStart','qEnd',
        'tName','tSize','tStart','tEnd',
        'blockCount','blockSizes','qStarts','tStarts']
    converters = []
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        super(Chain, self).__init__(*args,**kw)
    
    def convert(self):
        for key in self.__dict__:
            try:
                self.__dict__[key] = int(self.__dict__[key])
            except:
                pass
        
        if len(self.strand)==2:
            self.strand = self.strand[1]
        for key in ['blockSizes','qStarts','tStarts']:
            self.__dict__[key] = [int(x) for x in self.__dict__[key].split(',')[0:-1]]
    
    def projectOntoString(self, seq):
        extracted = []
        for tStart,blockSize in zip(self.tStarts, self.blockSizes):
            tEnd = tStart+3*blockSize-1
            if self.strand=='+':
                s = seq[tStart:tEnd+1]
            else:
                s = sequence.reverseComplement(seq)[tStart:tEnd+1]
            extracted.append(s)
        return extracted
    
    
def load_iter(iFileHandle, format='psl', **kw):
    """Return an iterator to the BLAT file.
    
    @param iFileHandle: Input filename or file.
    @param format: BLAT output format (optional; default: 'psl')
    """
    if not format in ['psl']:
        raise 'Only psl is currently supported.'
    
    iFile = smartopen(iFileHandle)
    
    skip = kw.pop('skip', 5)
    for i in xrange(skip):
        junk = iFile.readline()
    
    for line in iFile:
        if line:
            tokens = line.strip().split('\t')
            yield Chain(tokens, **kw)


def load(iFileHandle, format='psl', **kw):
    """Load BLAT alignment data.
    
    @param iFileHandle: Input filename or file.
    @param format: BLAT output format (only psl is currently supported)
    """
    data = []
    for chain in load_iter(iFileHandle, **kw):
        data.append(chain)
    return data
    
