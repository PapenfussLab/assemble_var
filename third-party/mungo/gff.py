"""
gff module
"""

import sys, re
import fasta
from mungoCore import *
from useful import smartopen, raiseDeprecated


groupSplitterPattern = re.compile('\s?;\s?')
kvPairSplitterPattern = re.compile('\s?=\s?|\s*')

def groupSplitter(group):
    return [x.strip() for x in groupSplitterPattern.split(group)]

def kvPairSplitter(kvPair):
    return kvPairSplitterPattern.split(kvPair)


def GffFile(iFileHandle, simple=False, **kw):
    """Factory function for Reader and Writer classes
    
    @param iFileHandle: GFF file name or object
    """
    if simple:
        return SimpleGffReader(iFileHandle, **kw)
    else:
        return GffReader(iFileHandle, **kw)


class SimpleGffReader(AbstractDataReader):
    def __init__(self, iFileHandle, **kw):
        """Constructor
        
        @param iFileHandle: Input file or name
        """
        super(SimpleGffReader, self).__init__(iFileHandle)
    
    def _generator(self):
        for line in self.iFile:
            line = line.strip()
            if line and line[0]!="#":
                tokens = line.split('\t')
                yield Feature(tokens)


class GffReader(AbstractDataReader):
    def __init__(self, iFileHandle, **kw):
        """Constructor
        
        @param iFileHandle: Input file or name
        """
        super(GffReader, self).__init__(iFileHandle)
    
    def _generator(self):
        lastName = None
        gene = []
        for line in self.iFile:
            line = line.strip()
            if line and line[0]!="#":
                tokens = line.split('\t')
                f = Feature(tokens)
                if f.name!=lastName:
                    if lastName!=None:
                        yield gene
                    gene = [f]
                    lastName = copy.copy(f.name)
                else:
                    gene.append(f)
                yield gene


class Feature(AbstractFeature):
    """GFF feature class
    
    @see: http://www.sanger.ac.uk/Software/formats/GFF/
    
    g = Feature(
        reference=,
        source=,
        type=,
        start=,
        end=,
        score=,
        strand=,
        phase=,
        group=
    )
    
    """
    
    attributes = ['reference','source','type','start','end','score',
        'strand','phase','group']
    converters = zip(['start','end','score'],[int,int,int])
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        """Constructor:
        
        @param args: GFF values
        @type args: list or dict
        @keyword reference: Reference
        @keyword source: Source
        @keyword type: GFF type (eg. mRNA, CDS, ...)
        @keyword start: Start
        @keyword end: End
        @keyword score: Score
        @keyword strand: Strand
        @keyword phase: Phase
        @keyword group: Group and other semi-colon separated key-value pairs
        """
        super(Feature, self).__init__(*args, **kw)
    
    def convert(self):
        super(Feature, self).convert()
        self.groupDict = self._groupToDict()
        self.name = self._extractName()
    
    def __repr__(self):
        d = {}
        for k,v in self.__dict__.iteritems():
            if v==None:
                v = '.'
            d[k] = v
        return self.format % d
    
    def _extractName(self):
        try:
            groupFields = groupSplitter(self.group)
            bits = kvPairSplitter(groupFields[0])
            name = ' '.join(bits[1:])
        except:
            name = ''
        return name
    
    def _groupToDict(self):
        d = {}
        try:
            groupKVs = groupSplitter(self.group)
            for kv in groupKVs:
                bits = kvPairSplitter(kv)
                k,v = bits[0], ' '.join(bits[1:])
                d[k] = v
        except:
            pass
        return d
    
    def groupToList(self):
        d = []
        try:
            groupKVs = groupSplitter(self.group)
            for kv in groupKVs:
                bits = kvPairSplitter(kv)
                k,v = bits[0], ' '.join(bits[1:])
                d.append([k,v])
        except:
            pass
        return d
    
    def groupFromList(self, kv):
        self.group = ' ; '.join(['%s %s' % (k,v) for k,v in kv])
    
    def extractName(self):
        raiseDeprecated('Feature.extractName(): Use name attribute')
        return self.name
    
    def groupToDict(self):
        raiseDeprecated('Feature.groupToDict(): Use groupDict attribute')
        return self.groupDict
    
    def extractGroupData(self):
        raiseDeprecated('Feature.extractGroupData(): Use groupDict attribute')
        return self.groupDict
    
    def getSequence(self, blastDb, padFivePrime=0, padThreePrime=0):
        start = max(1,self.start-padFivePrime)
        end = self.end+padThreePrime
        h,s = fasta.getSequence(blastDb, self.reference, start, end, self.strand)
        return h,s
    
    @staticmethod
    def fromDomain(d, refTransform=lambda d: 'chr%s' % d.accession, 
      source='Domains', _type='motif', aggregator='Motif', note=''):
        """Create a gff object from a hmmer domain (like) object
        
        Group field is constructed as:
        <aggregator> <domain.domain>; Note "<note>"; Evalue <d.eValue>
        
        @param d: hmmer domain object
        @param refTransform: reference transform function 
            (default: lambda d: 'chr%s' % d.accession)
        @param source: gff source value (default: 'Domains')
        @param _type: gff type value (default: 'motif')
        @param aggregator: aggregator (default: 'Motif')
        @param note: Note (default: '')
        """
        g = Feature(doConvert=False)
        g.reference = refTransform(d)
        g.source = source
        g.type = _type
        g.start = int(d.sStart)
        g.end = int(d.sEnd)
        g.strand = d.strand
        g.score = d.score
        group = ['%s %s' % (aggregator, d.domain)]
        if note:
            group.append('Note "%s"' % note)
        group.append('E-value %g' % float(d.eValue))
        g.group = '; '.join(group)
        return g
    
    @staticmethod
    def fromHSP(d, refTransform=lambda d: d.subjectId,
      source='blast', _type='HSP', aggregator='Match', note='',
      **kw):
        """Create a gff object from an HSP object"""
        
        g = Feature(doConvert=False)
        g.reference = refTransform(d)
        g.source = source
        g.type = _type
        g.start = d.sStart
        g.end = d.sEnd
        g.strand = d.strand()
        g.score = d.bitScore
        group = ['%s %s' % (aggregator, d.queryId)]
        if note: group.append('Note "%s"' % note)
        if kw:
            for k,v in kw.iteritems():
                group.append('%s %s' % (k,v))
        group.append('E-value %g' % d.eValue)
        g.group = '; '.join(group)
        return g
    
    @staticmethod
    def fromChain(chain, refTransform=lambda d: d.subjectId,
      source='blast', _type='match', aggregator='Match', note='',
      **kw):
        """Create a gff object from a chain of HSP objects"""
        
        min_eValue = 999
        cumScore = 0.0
        extrema = []
        for hsp in chain:
            extrema.append(hsp.sStart)
            extrema.append(hsp.sEnd)
            cumScore += hsp.bitScore
            min_eValue = min(min_eValue, hsp.eValue)
        
        hsp = chain[0]
        g = Feature(doConvert=False)
        g.reference = refTransform(hsp)
        g.source = source
        g.type = _type
        g.start = min(extrema)
        g.end = max(extrema)
        g.strand = hsp.strand()
        g.score = cumScore
        group = ['%s %s' % (aggregator, hsp.queryId)]
        if note: group.append('Note "%s"' % note)
        if kw:
            for k,v in kw.iteritems():
                group.append('%s %s' % (k,v))
        group.append('Min_E-value %g' % min_eValue)
        g.group = '; '.join(group)
        return g
