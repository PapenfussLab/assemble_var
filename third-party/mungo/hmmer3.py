"""
hmmer3 module
"""

from mungo.mungoCore import *
from mungo.useful import smartopen, extractRootName
import sys, re, warnings

hmmer2frame = {0: 1, 1: 2, 2: 3, 3: -1, 4: -2, 5: -3}
frame2hmmer = dict([(v,k) for k,v in hmmer2frame.iteritems()])


class Domain(AbstractFeature):
    """Domain feature class"""
    
    attributes = ["targetName","accession","targetLen","queryName",
        "accession","qlen","evalue","seqScore","seqBias","i","N",
        "c_evalue","i_evalue","domainScore","domainBias",
        "hmmStart","hmmEnd","alnStart","alnEnd", "envStart","envEnd",
        "acc","description"]
    converters = zip(
        ["hmmStart","hmmEnd","alnStart","alnEnd","envStart","envEnd",
            "seqScore","domainScore","evalue","c_evalue","i_evalue",
            "seqScore", "seqBias", "domainScore", "domainBias"],
        [int,int,int,int,int,int,float,float,float,float,float,float,float,float,float])
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        """Constructor"""
        super(Domain, self).__init__(*args, **kw)
        self.genomic = False
    
    def __repr__(self):
        d = {}
        for k,v in self.__dict__.iteritems():
            d[k] = v
        return self.format % d
    
    def getSequence(self, blastdb, getAll=False, convertAccession=lambda x: x):
        if getAll:
            start = 0
            end = 0
        else:
            start = self.alnStart
            end = self.alnEnd
        accession = convertAccession(self.accession)
        h,s = blast.getSequence(blastdb, accession, start, end)
        return h,s


class BlockSixFrameDomain(Domain):
    def toGenomic(self):
        """Convert from 6 frame to genomic coordinates."""
        prog = re.compile('\.|-|\:')
        chrom,blockStart,blockEnd,hmmerFrame = prog.split(self.targetName)
        blockStart = int(blockStart)
        blockEnd = int(blockEnd)
        L = blockEnd-blockStart+1
        hmmerFrame = int(hmmerFrame)
        frame = hmmer2frame[hmmerFrame]
        
        if frame>0:
            strand = '+'
        else:
            strand = '-'
        
        gStart,gEnd = convertSixFrameToGenomic(self.alnStart, self.alnEnd, frame, L)
        self.genomic = True
        self.targetName = chrom
        self.alnStart = gStart
        self.alnEnd = gEnd
        self.strand = strand


def convertSixFrameToGenomic(start, end, frame, L):
    """Convert 6 frame coords to genomic.
    
    @param start: Amino acid start coord
    @param end: Amino acid end coord
    @param frame: Frame
    @param L: Nucleotide seq length
    @return: (gStart, gEnd, strand)
    """
    if frame>=0:
        gStart = 3*(start-1)+(frame-1)+1
        gEnd = 3*(end-1)+(frame-1)+3
    else:
        gStart = L-(3*(start-1)+abs(frame)-1)
        gEnd = L-(3*(end-1)+abs(frame)+1)
    return gStart,gEnd


def HmmerFile(iFileHandle, **kw):
    "Factory function returning a HmmerFileReader"
    return DomainHitsReader(iFileHandle, **kw)


class DomainHitsReader(AbstractDataReader):
    def __init__(self, iFileHandle, seqType=None, eValueCutoff=None, scoreCutoff=None):
        self.seqType = seqType
        self.eValueCutoff = eValueCutoff
        self.scoreCutoff = scoreCutoff
        super(DomainHitsReader, self).__init__(iFileHandle)
    
    def _generator(self):
        """Return an iterator to a HMMer file."""
        if self.seqType in [Domain, BlockSixFrameDomain]:
            _Domain = self.seqType
        elif self.seqType=='SixFrame':
            _Domain = SixFrameDomain
        elif self.seqType=='BlockSixFrame':
            _Domain = BlockSixFrameDomain
        else:
            _Domain = Domain
        
        for line in self.iFile:
            line = line.strip()
            if line[0]=="#": continue
            tokens = line.split()
            d = _Domain(dict(zip(Domain.attributes, tokens)))
            if (self.eValueCutoff and d.eValue>self.eValueCutoff) or \
              (self.scoreCutoff and d.score<self.scoreCutoff): continue
            yield d

