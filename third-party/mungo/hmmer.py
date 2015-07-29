"""
hmmer module
"""

from __future__ import print_function

from mungo.mungoCore import *
import blast, sequence
from mungo.useful import smartopen, extractRootName, ClassFromDict, warnDeprecated
import sys, re, warnings


hmmer2frame = {0: 1, 1: 2, 2: 3, 3: -1, 4: -2, 5: -3}
frame2hmmer = dict([(v,k) for k,v in hmmer2frame.iteritems()])


def HmmerFile(iFileHandle, **kw):
    "Factory function returning a HmmerFileReader"
    return HmmerReader(iFileHandle, **kw)


class HmmerReader(AbstractDataReader):
    def __init__(self, iFileHandle, seqType=None, eValueCutoff=None, scoreCutoff=None):
        super(HmmerReader, self).__init__(iFileHandle)
        self.seqType = seqType
        self.eValueCutoff = eValueCutoff
        self.scoreCutoff = scoreCutoff
    
    def _generator(self):
        """Return an iterator to a HMMer file."""
        if self.seqType in [Domain, SixFrameDomain, BlockSixFrameDomain, OrfDomain, OrfDomain2]:
            _Domain = self.seqType
        elif self.seqType=='SixFrame':
            _Domain = SixFrameDomain
        elif self.seqType=='BlockSixFrame':
            _Domain = BlockSixFrameDomain
        elif self.seqType=='ORFs':
            _Domain = OrfDomain
        else:
            _Domain = Domain
        
        startToken = '^Parsed for domains'
        endToken = '^Alignments of top-scoring domains'
        abortToken = '\[no hits above thresholds\]'
        
        startRegex = re.compile(startToken)
        if not jumpToMatch(self.iFile, startRegex):
            raise Exception('No match found. File may be empty.')
        
        # 3. Parse domain details
        line = self.iFile.next()
        line = self.iFile.next()
        endRegex = re.compile(endToken)
        abortRegex = re.compile(abortToken)
        domains = []
        for line in self.iFile:
            line = line.strip()
            if endRegex.match(line) or abortRegex.match(line):
                break
            elif not line:
                continue
            tokens = line.split()
            d = _Domain(dict(zip(Domain.attributes[1:], tokens)))
            if (self.eValueCutoff and d.eValue>self.eValueCutoff) or \
              (self.scoreCutoff and d.score<self.scoreCutoff): continue
            yield d


class PfamReader(AbstractDataReader):
    def __init__(self, iFileHandle, eValueCutoff=None, scoreCutoff=None):
        super(PfamReader, self).__init__(iFileHandle)
        self.eValueCutoff = eValueCutoff
        self.scoreCutoff = scoreCutoff
    
    def _generator(self):
        pass


class Domain(AbstractFeature):
    """Domain feature class"""
    
    attributes = ['domain', 'accession', 'count', 'sStart', 'sEnd', 
        'sCode', 'qStart', 'qEnd', 'qCode', 'score', 'eValue']
    converters = zip(
        ['qStart','qEnd','sStart','sEnd','score','eValue'],
        [int,int,int,int,float,float])
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        """Constructor:
        
        @param args: HMMer field values
        @type args: list, dict, Domain
    
        Optional keywords:
        @keyword domain: Domain name
        @keyword accession: Subject name
        @keyword count: Id/total hits on subject
        @keyword sStart:
        @keyword sEnd:
        @keyword sCode:
        @keyword qStart:
        @keyword qEnd:
        @keyword qCode:
        @keyword score: Bit score
        @keyword eValue:
        """
        super(Domain, self).__init__(*args, **kw)
        self.genomic = False
    
    def __repr__(self):
        d = {}
        for k,v in self.__dict__.iteritems():
            d[k] = v
        return self.format % d
    
    def getTokens(self):
        return [self.__dict__[key] for key in self.attributes]
    
    def addAttribute(self, attribute, default=None):
        self.attributes.append(attribute)
        self.format = self.format + '\t%%(%s)s' % attribute
        self.__dict__[attribute] = default
    
    def addStrandAttribute(self, strand=None):
        self.addAttribute('strand', strand)
    
    def swapStartEnd(self):
        if self.sStart>self.sEnd:
            self.sStart,self.sEnd = self.sEnd,self.sStart
    
    def getSequence(self, blastdb, getAll=False, convertAccession=lambda x: x):
        if getAll:
            start = 0
            end = 0
        else:
            start = self.sStart
            end = self.sEnd
        accession = convertAccession(self.accession)
        h,s = blast.getSequence(blastdb, accession, start, end)
        return h,s
    
    @staticmethod
    def fromGenomic(tokens):
        strand = tokens[-1]
        d = Domain(tokens[0:-1])
        d.genomic = True
        d.addStrandAttribute(strand)
        return d


class OrfDomain(Domain):
    def toGenomic(self, orfStart, orfStrand, doSwapStartEnd=True):
        """Convert from ORF to genomic coordinates."""
        self.genomic = True
        self.sStart,self.sEnd = convertOrfToGenomic(
            self.sStart, self.sEnd, orfStrand, orfStart)
        self.addStrandAttribute(orfStrand)
        if doSwapStartEnd:
            self.swapStartEnd()


class OrfDomain2(Domain):
    "ORF domain class for use with my ORF files"
    def toGenomic(self, doSwapStartEnd=True):
        """Convert from ORF to genomic coordinates."""
        self.genomic = True
        o = parseOrfHeader(self.accession)
        self.sStart,self.sEnd = convertOrfToGenomic(
            self.sStart, self.sEnd, o.strand, o.start)
        self.addStrandAttribute(o.strand)
        if doSwapStartEnd:
            self.swapStartEnd()


class SixFrameDomain(Domain):
    def toGenomic(self, L, doSwapStartEnd=True):
        """Convert from 6 frame to genomic coordinates.
        
        @param L: Length of DNA sequence.
        """
        self.genomic = True
        o = parseSixFrameHeader(self.accession)
        self.sStart,self.sEnd = convertSixFrameToGenomic(
            self.sStart, self.sEnd, o.frame, L)
        self.accession = o.name
        self.strand = o.strand
        self.addStrandAttribute(o.strand)
        if doSwapStartEnd:
            self.swapStartEnd()
    
    def toBlockCoords(self, L=1e99, blockSize=5000000, delimiter='.'):
        self.accession, self.sStart, self.sEnd = \
            blast.genomeToBlock(
                self.accession, self.sStart, self.sEnd, L=L, 
                blockSize=blockSize, delimiter=delimiter)
    
    def getSequenceFromString(self, seq):
        s = seq[self.sStart-1:self.sEnd]
        if self.strand=='-':
            s = sequence.reverseComplement(s)
        return s
    
    def getSequence(self, blastDb, padFivePrime=0, padThreePrime=0):
        if self.genomic:
            start = max(1,self.sStart-padFivePrime)
            end = self.sEnd+padThreePrime
            h,s = blast.getSequence(blastDb, self.accession, start, end, self.strand)
        else:
            raise Exception('You must call the toGenomic method first.')
        return h,s


class BlockSixFrameDomain(Domain):
    def toGenomic(self, relative=False, doSwapStartEnd=True, relDelimiter=':'):
        """Convert from 6 frame to genomic coordinates."""
        self.genomic = True
        chrom,blockStart,blockEnd,gStart,gEnd,strand = \
            convertBlockSixFrameToGenomic(
            self.accession, self.sStart, self.sEnd)
        
        if relative:
            self.accession = '%s%s%i-%i' % (chrom,relDelimiter,blockStart,blockEnd)
            self.sStart = gStart
            self.sEnd = gEnd
        else:
            self.accession = chrom
            self.sStart = blockStart + gStart - 1
            self.sEnd = blockStart + gEnd - 1
        self.addStrandAttribute(strand)
        if doSwapStartEnd:
            self.swapStartEnd()


class GenomicDomain(AbstractFeature):
    """GenomicDomain feature class"""
    
    attributes = ['domain', 'accession', 'count', 'sStart', 'sEnd', 
        'sCode', 'qStart', 'qEnd', 'qCode', 'score', 'eValue']
    converters = zip(
        ['qStart','qEnd','sStart','sEnd','score','eValue'],
        [int,int,int,int,float,float])
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        """Constructor:
        
        @param args: HMMer field values
        @type args: list, dict, Domain
    
        Optional keywords:
        @keyword domain: Domain name
        @keyword accession: Subject name
        @keyword count: Id/total hits on subject
        @keyword sStart:
        @keyword sEnd:
        @keyword sCode:
        @keyword qStart:
        @keyword qEnd:
        @keyword qCode:
        @keyword score: Bit score
        @keyword eValue:
        @keyword strand:
        """
        
        super(GenomicDomain, self).__init__(*args, **kw)
    
    def toDict(self):
        return self.__dict__
    
    def toList(self):
        return self.__dict__.items()
    
    def __repr__(self):
        try:
            d = {}
            for k,v in self.__dict__.iteritems():
                d[k] = v
            return self.format % d
        except:
            return str(self.__dict__)
    
    def toBlockCoords(self, L=1e99, blockSize=5000000, delimiter='.'):
        self.accession, self.sStart, self.sEnd = \
            blast.genomeToBlock(
                self.accession, self.sStart, self.sEnd, L=L, 
                blockSize=blockSize, delimiter=delimiter)
    
    def getSequence(self, blastDb, padFivePrime=0, padThreePrime=0):
        start = max(1,self.sStart-padFivePrime)
        end = self.sEnd+padThreePrime
        h,s = blast.getSequence(blastDb, self.accession, start, end, self.strand)
        return h,s


def loadGenomicDomains(filename):
    data = []
    gene = []
    for line in open(filename):
        line = line.strip()
        if not line:
            continue
        elif line[0] in ['#', '>']:
            if gene:
                data.append(gene)
            gene = []
        else:
            tokens = line.split('\t')
            d = GenomicDomain(tokens)
            gene.append(d)
    data.append(gene)
    return data


def jumpToMatch(iFile, regex):
    """Jump to regex match in file.
    
    @param iFile: File object
    @param regex: Compiled regex object
    @return: True if successful, False otherwise
    """
    for line in iFile:
        if regex.match(line):
            return True
    return False


def extractUptoMatch(iFile, regex):
    """Extract up to regex match from file.
    
    @param iFile: File object
    @param regex: Compiled regex object
    @return: string
    """
    block = []
    for line in iFile:
        if regex.match(line):
            break
        else:
            block.append(line.rstrip())
    return block


def parseSixFrameHeader(header):
    """Parse a 6 frame header (from translate or python).
    
    @param header: Six frame header "<name>:<frame>" or "<name>.<start>-<end>:<frame>"
                   (assumes input frame is hmmer frame (0-5)).
    @return: a simple class with attributes name, start, end, strand and frame.
    """
    
    header = header.strip()
    regex = re.compile(
        '(?P<name>\w+)([\.|:](?P<start>\d+)[-|,](?P<end>\d+))?:(?P<frame>[0-5])')
    rs = regex.search(header)
    d = rs.groupdict()
    
    d['frame'] = hmmer2frame[int(d['frame'])]
    
    if d['frame']>0:
        d['strand'] = '+'
    else:
        d['strand'] = '-'
    
    try:
        d['start'] = int(d['start'])
        d['end'] = int(d['end'])
    except:
        pass
    
    return ClassFromDict(d)


def parseOrfHeader(header):
    """Parse an ORF header (from extractORFs.py).
    
    @param header: ORF header "<name>.<orfId>.<start>-<end>  Length=<length>"
                   (Length optional).
    @return: a simple class with attributes name, start, end, strand and length.
    """
    regex = re.compile(
        '(?P<name>\w+)\.(?P<orfId>\d+)\.(?P<start>\d+)-(?P<end>\d+)(\SLength=(?P<length>\d+))?')
    rs = regex.match(header.strip())
    d = rs.groupdict()
    try:
        d['start'] = int(d['start'])
        d['end'] = int(d['end'])
        d['length'] = int(d['length'])
    except:
        pass
    
    if d['start']>d['end']:
        d['strand'] = '-'
    else:
        d['strand'] = '+'
    
    return ClassFromDict(d)


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


def convertBlockSixFrameToGenomic(block, start, end):
    """Convenience function that takes block 6 frame coords 
    (block,start,end), extracts the block start/end and frame
    and converts them to genomic coords
    
    ie. 
    
    chrom.blockStart-blockEnd:frame aaStart aaEnd or
    chrom:blockStart-blockEnd:frame aaStart aaEnd
    --> chrom,blockStart,blockEnd,gStart,gEnd,strand
    
    @param block: Block accession ("<name>.<blockStart>-<blockEnd>:<frame>")
    @param start: Domain start
    @param end: Domain end
    @return: (chrom, blockStart, blockEnd, gStart, gEnd, strand)
    """
    prog = re.compile('\.|-|\:')
    tokens = prog.split(block)
    if len(tokens)==2:
        chrom,hmmerFrame = tokens
        blockStart = 0
        blockEnd = 0
        L = 0
    elif len(tokens)==4:
        chrom,blockStart,blockEnd,hmmerFrame = tokens
        blockStart = int(blockStart)
        blockEnd = int(blockEnd)
        L = blockEnd-blockStart+1
    else:
        print(tokens, file=sys.stderr)
        raise Exception("Don't know what to do")
    
    hmmerFrame = int(hmmerFrame)
    frame = hmmer2frame[hmmerFrame]
    if frame>0:
        strand = '+'
    else:
        strand = '-'
    gStart,gEnd = convertSixFrameToGenomic(start, end, frame, L)
    return chrom,blockStart,blockEnd,gStart,gEnd,strand


def convertGenomicToBlockCoords(domain, chrLen, blockSize=5000000, delimiter='.'):
    domain.accession, domain.sStart, domain.sEnd = \
        blast.genomeToBlock(
            domain.accession, domain.sStart, domain.sEnd, 
            L=chrLen, blockSize=blockSize, delimiter=delimiter)
    return domain


def convertOrfToGenomic(start, end, strand, orfStart):
    """Convert domain coordinates in ORF to genomic.
    
    @param start: Domain start coord
    @param end: Domain end coord
    @param strand: Strand
    @param orfStart: ORF start coord
    @return: (gStart, gEnd)
    """
    if strand=='+':
        gStart = orfStart + 3*(start-1)
        gEnd = orfStart + 3*(end-1) + 2
    else:
        gStart = orfStart - 3*(start-1)
        gEnd = orfStart - 3*(end-1) - 2
    return gStart, gEnd


def loadDomains(iFileHandle):
    """Load hmmer domain results.
    
    @param iFileHandle: Input file or filename
    @param seqType: Type of sequence searched 
        [None (default), 'SixFrame', 'BlockSixFrame' or 'ORFs']
    @param eValueCutoff: E-value threshold (default None)
    @param scoreCutoff: Score threshold (default None)
    @return: list of domains
    """
    
    domains = []
    for d in HmmerFile(iFileHandle):
        domains.append(d)
    return domains
