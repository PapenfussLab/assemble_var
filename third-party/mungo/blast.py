"""
blast module

Requirements:
    - getSequence uses BLAST's fastacmd
"""

import os, sys, copy, re
from mungoCore import *
import fasta, sequence
from useful import smartopen


def sign_to_strand(x):
    if x<0:
        return '-'
    else:
        return '+'

def test_to_strand(test):
    if test:
        return '+'
    else: 
        return '-'


class HSP(AbstractFeature):
    """HSP feature"""
    
    attributes = ['queryId','subjectId','pcId','alignLength',
        'matches','mismatches','qStart','qEnd','sStart','sEnd','eValue',
        'bitScore']
    converters = [
        ('pcId', float),
        ('alignLength', int),
        ('matches', int),
        ('mismatches', int),
        ('qStart', int),
        ('qEnd', int),
        ('sStart', int),
        ('sEnd', int),
        ('eValue', float),
        ('bitScore', float)
    ]
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        super(HSP,self).__init__(*args,**kw)
        if self.sStart>self.sEnd:
            self.sStart,self.sEnd = self.sEnd,self.sStart
            self._strand = '-'
        else:
            self._strand = '+'
    
    def __repr__(self):
        self.forcePlusStrand()
        output = self.format % self.__dict__
        self.forceCorrectStrand()
        return output
    
    def __eq__(self,x):
        if self.__dict__==x.__dict__:
            return True
        else:
            return False
    
    def __hash__(self):
        return hash("%(queryId)s %(subjectId)s:%(sStart)i-%(sEnd)i" % self.__dict__)
    
    def addField(self, field, default=None):
        if not field in self.attributes:
            self.attributes.append(field)
            self.format = self.format + '\t%%(%s)s' % field
            self.__dict__[field] = default
    
    def strand(self):
        strand = test_to_strand(
            self._strand==sign_to_strand(self.sEnd-self.sStart))
        return strand
    
    def swapStartEnd(self):
        self.sStart,self.sEnd = self.sEnd,self.sStart
        self._strand = {'+': '-', '-': '+'}[self._strand]
    
    def forcePlusStrand(self):
        if self._strand=='-':
            self.sStart,self.sEnd = self.sEnd,self.sStart
            self._strand = '+'
    
    def forceCorrectStrand(self):
        if self._strand=='+' and self.sStart<self.sEnd:
            self.sStart,self.sEnd = self.sEnd,self.sStart
            self._strand = '-'
    
    def convertBlockToGenomeCoords(self):
        """Convert block to genome coordinates."""
        try:
            for delimiter in [':', '.', '-']:
                self.subjectId = self.subjectId.replace(delimiter, ' ')
            tokens = self.subjectId.split()
            self.subjectId = tokens[0]
            self.sStart += int(tokens[1])-1
            self.sEnd += int(tokens[1])-1
        except Exception, e:
            print e
    
    def convertGenomeToBlockCoords(self, L, blockSize=5000000):
        """Convert genome to block coordinates.
        
        @param L: Length of chromosome
        @param blockSize: Length of block (Default = 5000000)
        
        """
        self.subjectId,self.sStart,self.sEnd = fasta.genomeToBlock(
            self.subjectId, self.sStart, self.sEnd, L, blockSize=blockSize)
    
    def getSequence(self, blastDb, accessionConverter=lambda x: x, 
      padFivePrime=0, padThreePrime=0, getAll=True):
        if getAll:
            start = 0
            end = 0
        else:
            start = max(1,self.sStart-padFivePrime)
            end = self.sEnd+padThreePrime
        accession = accessionConverter(self.subjectId)
        h,s = getSequence(blastDb, accession, start, end, self.strand())
        return h,s


def BlastFile(iFileHandle, multi=False, **kw):
    """Factory function for Reader and Writer classes
    
    @param iFileHandle: BLAST file name or object
    """
    if multi:
        return MultiBlastReader(iFileHandle, **kw)
    else:
        return BlastReader(iFileHandle, **kw)


class BlastReader(AbstractDataReader):
    def __init__(self, iFileHandle, eValueCutoff=None, **kw):
        """Constructor
        
        @param iFileHandle: Input file or name
        """
        super(BlastReader, self).__init__(iFileHandle)
        self.eValueCutoff = eValueCutoff
    
    def iterBestHits(self):
        oldQueryId = None
        for hsp in self:
            if hsp.queryId!=oldQueryId:
                yield hsp
                oldQueryId = copy.copy(hsp.queryId)
    
    def bestHits(self):
        data = []
        for bh in self.iterBestHits():
            data.append(bh)
        return bh
    
    def firstHit(self):
        for hsp in self:
            return hsp
    
    def _generator(self):
        for line in self.iFile:
            line = line.strip()
            if line and line[0]!='#':
                tokens = line.split('\t')
                hsp = HSP(tokens)
                if not self.eValueCutoff or hsp.eValue<=self.eValueCutoff:
                    yield HSP(tokens)


class MultiBlastReader(BlastReader):
    """
    Subclass of BlastReader for parsing the results of aligning
    a multi-fasta file to a blast database.
    """
    
    def bestHits(self):
        hsps = []
        oldQueryId = None
        for hsp in self:
            if hsp.queryId!=oldQueryId and not oldQueryId is None:
                if len(hsps)>1:
                    yield hsps[0], hsps[1], len(hsps)
                else:
                    yield hsps[0], None, len(hsps)
                hsps = []
            hsps.append(hsp)
            oldQueryId = copy.copy(hsp.queryId)
        
        if len(hsps)>1:
            yield hsps[0], hsps[1], len(hsps)
        elif len(hsps)==1:
            yield hsps[0], None, len(hsps)
    
    def groupByQuerySeq(self):
        hsps = []
        oldQueryId = None
        for hsp in self:
            if hsp.queryId!=oldQueryId and not oldQueryId is None:
                yield hsps
                hsps = []
            hsps.append(hsp)
            oldQueryId = copy.copy(hsp.queryId)
        yield hsps


class NotFoundException(Exception):
    pass


def getSequence(blastDb, accession, start=0, end=0, strand='+', padding=0, debug=False):
    """Load a sequence from a BLAST database.
    
    @param blastDb: BLAST database
    @param accession: Accession name
    @param start: Start coordinate (Default: 0, extract from start of sequence)
    @param end: End coordinate (Default: 0, extract to the end of sequence)
    @param strand: Strand (Default: '+')
    @param padding: Sequence padding (Default: 0)
    @returns: (header,seq)
    """
    if start>end: start,end = end,start
    
    cmd = 'fastacmd -d %s -s "%s" -L %i,%i' % (blastDb,accession,start,end)
    if debug: print cmd
    p = os.popen(cmd)
    header = p.readline()[1:].strip()
    if not header:
        raise Exception('BLAST failure')
    seq = []
    for line in p:
        seq.append(line.strip())
    seq = ''.join(seq)
    
    if not seq:
        print blastDb, accession
        raise NotFoundException()
    
    if strand=='-':
        seq = sequence.reverseComplement(seq)
    
    return header,seq


class TooManyBlocks(Exception):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        print chrom,start,end


def genomeToBlock(chrom, start, end, L=1000000000, blockSize=5000000, delimiter='.'):
    """Transform genomic to block coordinates.
    
    @param chrom: Chromosome
    @param start: Start coordinate
    @param end: End coorinate
    @param L: Chromosome length
    @param blockSize: Block size (Default: 5000000)
    @returns: (block, relStart, relEnd)
    """
    strand = '+'
    if start>end:
        strand = '-'
        start,end = end,start
    
    blockStart = 1 + blockSize*int(start/blockSize)
    blockEnd = min(L, blockStart+blockSize-1)
    block = '%s%s%i-%i' % (chrom, delimiter, blockStart, blockEnd)
    
    relStart = start % blockSize
    relEnd = end % blockSize
    if strand=='-':
        relStart,relEnd = relEnd,relStart
    
    return block,relStart,relEnd


def genomeToBlocks(chrom, start, end, L, blockSize=5000000):
    """Transform genomic to block coordinates 
    (chrom.blockStart-blockEnd:start-end).
    
    @param chrom: Chromosome
    @param start: Start coordinate
    @param end: End coorinate
    @param L: Chromosome length
    @param blockSize: Block size (Default: 5000000)
    @returns: a tuple of lists (blocks,startInBlocks,endInBlocks).
    """
    strand = '+'
    if start>end:
        strand = '-'
        start,end = end,start
    
    startBlock = 1 + blockSize*int(start/blockSize)
    endBlock = min(L, blockSize*(1+int(end/blockSize)))
    
    blocks = []
    i = 0 ; imax = 200
    s = startBlock
    while s<endBlock:
        e = min(L, s+blockSize-1)
        blocks.append('%s.%i-%i' % (chrom,s,e))
        s += blockSize
        i += 1
        if i>imax:
            raise TooManyBlocks(chrom,start,end)
    
    startInBlocks = [start % blockSize]
    endInBlocks = [end % blockSize]
    
    for i in xrange(len(blocks)-1):
        startInBlocks.append(0)
        endInBlocks.insert(0,0)
    
    if strand=='-':
        startInBlocks,endInBlocks = endInBlocks,startInBlocks
    return blocks,startInBlocks,endInBlocks


def blockToGenome(block, startInBlock, endInBlock):
    """Transform block to genomic coordinates.
    
    @param block: Block (chrom.blockStart-blockEnd)
    @param startInBlock: Start coordinate in block
    @param endInBlock: End coorinate in block
    @returns: (chrom, gStart, gEnd)
    """
    prog = re.compile('[\.|:]|-')
    chrom,blockStart,blockEnd = prog.split(block)
    gStart = int(blockStart)+startInBlock-1
    gEnd = int(blockStart)+endInBlock-1
    return chrom,gStart,gEnd
