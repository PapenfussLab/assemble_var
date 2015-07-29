"""
bowtie
"""

from mungoCore import AbstractDataReader
from useful import smartopen


class Feature(object):
    def __init__(self, tokens):
        self.tokens = tokens
        self.name = tokens[0]
        self.strand = tokens[1]
        self.chrom = tokens[2]
        self.start = int(tokens[3]) # 1-indexed
        # self.seq = tokens[4]
        # self.qual = tokens[5]
        # self.reserved = tokens[6]
        # self.decript = tokens[7]
    
    def __repr__(self):
        return "\t".join(self.tokens)


class PairedEnds(object):
    def __init__(self, feature1, feature2):
        self.ends = [feature1, feature2]
        self.sort()
    
    def sort(self):
        if self.ends[0].chrom==self.ends[1].chrom \
        and self.ends[0].start>self.ends[1].start:
            self.ends.reverse()
    
    def __getitem__(self, i):
        return self.ends[i]
    
    def fromTokens(tokens1, tokens2):
        self.ends = [Feature(tokens1), Feature(tokens2)]


class BowtieReader(AbstractDataReader):
    def __init__(self, iFileHandle, paired=False):
        super(BowtieReader, self).__init__(iFileHandle)
        self.paired = paired
    
    def _generator(self):
        if not self.paired:
            for line in self.iFile:
                tokens = line.strip().split('\t')
                yield Feature(tokens)
        elif self.paired:
            data = {}
            for line in self.iFile:
                tokens = line.strip().split('\t')
                f = Feature(tokens)
                name = f.name.split("/")[0]
                try:
                    data[name].append(f)
                    pe = data.pop(name)
                    yield PairedEnds(pe[0], pe[1])
                except KeyError:
                    data[name] = [f]
        else:
            raise Exception("BowtieReader: paired!=True/False")


class BowtieTwoFilePairedEndReader(object):
    def __init__(self, iFileHandle1, iFileHandle2):
        self.iFile1 = smartopen(iFileHandle1)
        self.iFile2 = smartopen(iFileHandle2)
    
    def __iter__(self):
        self._iter = self._generator()
        return self
    
    def next(self):
        for x in self._iter:
            return x
        raise StopIteration
    
    def _generator(self):
        data = {}
        go = True
        while go:
            for iFile in [self.iFile1, self.iFile2]:
                line = iFile.readline()
                if not line: 
                    go = False
                    continue
                
                tokens = line.strip().split("\t")
                f = Feature(tokens)
                name = f.name.split("/")[0]
                try:
                    data[name].append(f)
                    pe = data.pop(name)
                    yield PairedEnds(pe[0], pe[1])
                except KeyError:
                    data[name] = [f]
