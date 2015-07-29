"""
bed module
"""

from useful import smartopen
from mungoCore import *


class Feature(AbstractFeature):
    "BED feature"
    
    attributes = ['chrom','chromStart','chromEnd','name','score','strand',
        'thickStart','thickEnd','reserved','blockCount','blockSizes','blockStarts']
    requiredAttributes = ['chrom','chromStart','chromEnd']
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        """Feature constructor.
        
        @param args: Optional BED attribute values (same order as attributes)
        @type args: list
        
        Optional keywords:
        @keyword chrom: Chromosome
        @keyword chromStart: Start indexed from 0
        @keyword chromEnd: End indexed from 1
        @keyword name: Name
        @keyword score: Score
        @keyword strand: Strand
        @keyword thickStart: Eg. start of coding region
        @keyword thickEnd: Eg. end of coding region
        @keyword reserved:
        @keyword blockCount: Number of blocks
        @keyword blockSizes: Block lengths
        @keyword blockStarts: Start of block positions
        
        Only chrom, chromStart & chromEnd are required attributes 
        (either as elements of tokens or as keywords).
        """
        super(Feature, self).__init__(*args,**kw)
        self.checkConstraints()

    def checkConstraints(self):
        for attrib in Feature.requiredAttributes:
            if not self.__dict__[attrib]:
                raise '%s is a required attributes' % attrib
    
    def convert(self):
        for attrib in ['chromStart','chromEnd','thickStart','thickEnd','score']:
            try:
                self.__dict__[attrib] = int(self.__dict__[attrib])
            except:
                pass
        try:
            if self.blockStarts:
                self.blockStarts = [int(x) for x in self.blockStarts.split(',')[0:-1]]
        except:
            pass
        
        try:
            if self.blockSizes:
                self.blockSizes = [int(x) for x in self.blockSizes.split(',')[0:-1]]
        except:
            pass
    
    def __repr__(self):
        output = []
        for attribute in self.attributes:
            try:
                value = self.__dict__[attribute]
                output.append(str(value))
            except KeyError:
                pass
        return "\t".join(output)


def load(iFilename, offset=0):
    """Load a BED file.

    @param iFilename: Input filename or file.
    @param offset: Offset subtracted from positions (Default: 0).
    @return: List of features.
    """
    iFile = smartopen(iFilename)
    data = []
    for line in iFile:
        line = line.strip()
        if line and line[0]!='#':
            tokens = line.split('\t')
            f = Feature(tokens)
            
            try:
                f.chromStart -= offset
                f.chromEnd -= offset
                f.thickStart -= offset
                f.thickEnd -= offset
            except:
                pass
            
            data.append(f)
    return data


def load_iter(iFilename, offset=0):
    """Load a BED file.

    @param iFilename: Input filename or file.
    @param offset: Offset subtracted from positions (Default: 0).
    @return: List of features.
    """
    iFile = smartopen(iFilename)
    for line in iFile:
        line = line.strip()
        if line and line[0]!='#':
            tokens = line.split('\t')
            f = Feature(tokens)
            
            try:
                f.chromStart -= offset
                f.chromEnd -= offset
                f.thickStart -= offset
                f.thickEnd -= offset
            except:
                pass
            
            yield f

