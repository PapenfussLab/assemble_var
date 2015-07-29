"""
generic feature module
"""

from mungoCore import *
from useful import smartopen


def GenericFile(iFileHandle, **kw):
    return GenericReader(iFileHandle, **kw)


class GenericReader(AbstractDataReader):
    def __init__(self, iFileHandle, delimiter='\t', header=False, blockCol=None):
        super(GenericReader, self).__init__(iFileHandle)
        self.delimiter = delimiter
        self.header = header
        self.blockCol = blockCol
    
    def _generator(self):
        if self.header:
            junk = self.iFile.next()
        
        if self.blockCol:
            last = None
            block = []
            for line in self.iFile:
                tokens = line[0:-1].split(self.delimiter)
                if tokens[self.blockCol]!=last:
                    if last:
                        yield block
                    block = [tokens]
                    last = copy.copy(tokens[self.blockCol])
                else:
                    block.append(tokens)
            yield block
        else:
            for line in self.iFile:
                yield line.strip().split(self.delimiter)


class Feature(AbstractFeature):
    """Generic feature class"""
    
    def __init__(self, attributes, tokens, converters=None):
        self.attributes = attributes
        self.format = attributesToFormat(attributes)
        if converters:
            self.converters = converters
        else:
            self.converters = []
        super(Feature, self).__init__(tokens)


def load_iter(iFileHandle, skip=0, attributes=None, converters=None, debug=False):
    iFile = open(iFileHandle)
    for i in xrange(skip): iFile.next()
    if not attributes:
        attributes = iFile.readline().strip().split('\t')
    if debug: print attributes
    for line in iFile:
        if line:
            tokens = line.strip().split('\t')
            if debug: print tokens
            f = Feature(attributes, tokens, converters=converters)
            if debug: print f
            yield f

def load(iFileHandle, **kw):
    return [f for f in load_iter(iFileHandle, **kw)]
