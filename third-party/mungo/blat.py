"""
blat module
"""

from mungoCore import *
import sequence
from useful import smartopen, DSU


def string2list(string):
    return [int(x) for x in string.split(',')[:-1]]


def list2string(array):
    return ','.join([str(x) for x in array]) + ','


class Chain(AbstractFeature):
    """Blat chain feature"""
    
    attributes = [
        'match','mismatch','repMatch','Ns',
        'qGapCount','qGapBases',
        'tGapCount','tGapBases','strand',
        'qName','qSize','qStart','qEnd',
        'tName','tSize','tStart','tEnd',
        'blockCount','blockSizes','qStarts','tStarts']
    
    converters = [
        ('blockSizes', string2list),
        ('qStarts', string2list),
        ('tStarts', string2list)
    ]
    
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        super(Chain, self).__init__(*args,**kw)
    
    def __repr__(self):
        for key in ['blockSizes','qStarts','tStarts']:
            self.__dict__[key] = list2string(self.__dict__[key])
        out = self.format % self.__dict__
        for key in ['blockSizes','qStarts','tStarts']:
            self.__dict__[key] = string2list(self.__dict__[key])
        return out
    
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


def BlatFile(iFileHandle, multi=False, **kw):
    """Factory function for Reader and Writer classes
    
    @param iFileHandle: BLAT file name or object
    """
    if multi:
        return MultiBlatReader(iFileHandle, **kw)
    else:
        return BlatReader(iFileHandle, **kw)


class BlatReader(AbstractDataReader):
    def __init__(self, iFileHandle, eValueCutoff=None, **kw):
        """Constructor
        
        @param iFileHandle: Input file or name
        """
        super(BlatReader, self).__init__(iFileHandle)
    
    def _generator(self):
        for i in xrange(5):
            junk = self.iFile.next()
        
        for line in self.iFile:
            line = line.strip()
            tokens = line.split('\t')
            yield Chain(tokens)


class MultiBlatReader(BlatReader):
    def iterSorted(self):
        chains = []
        oldQName = None
        for chain in self:
            if chain.qName!=oldQName and not oldQName is None:
                chains.sort(key=lambda c: c.mismatch-c.match)
                yield chains
                chains = []
            chains.append(chain)
            oldQName = copy.copy(chain.qName)
        chains.sort(key=lambda c: c.mismatch-c.match)
        yield chains
    
    def iterBestHits(self):
        chains = []
        oldQName = None
        for chain in self:
            if chain.qName!=oldQName and not oldQName is None:
                chains.sort(key=lambda c: c.mismatch-c.match)
                yield chains[0]
                chains = []
            chains.append(chain)
            oldQName = copy.copy(chain.qName)
        chains.sort(key=lambda c: c.mismatch-c.match)
        yield chains[0]
