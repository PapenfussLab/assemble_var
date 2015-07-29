"""
TRANSFAC module
"""

from useful import smartopen
from mungoCore import *


class Match(AbstractFeature):
    "TRANSFAC match feature"
    
    attributes = ['matrix','start','strand','score','coreScore','seq','name']
    header = '\t'.join(attributes)
    converters = []
    format = attributesToFormat(attributes)
    
    def __init__(self,  *args, **kw):
        super(Match, self).__init__(*args,**kw)
    
    def convert(self):
        self.start = int(self.start)
        self.score = float(self.score)
        self.coreScore = float(self.coreScore)


def load(iFileHandle, skip=5, splitOn=None):
    """Load TRANSFAC match output.
    
    Arguments:
    iFileHandle -- Input file or filename.
    
    """
    data = []
    for match in load_iter(iFileHandle, skip, splitOn):
        data.append(match)
    return data


def load_iter(iFileHandle, skip=5, splitOn=None):
    """Load TRANSFAC match output.
    
    Arguments:
    iFileHandle -- Input file or filename.
    
    """
    iFile = smartopen(iFileHandle)
    for i in xrange(skip):
        iFile.next()
    
    for line in iFile:
        tokens = line.strip().split(splitOn)
        yield Match(tokens)
