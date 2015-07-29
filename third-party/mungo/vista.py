"""
vista2 module
"""

from mungoCore import *
from useful import smartopen


def convertStrandToPM(strand):
    if strand in ['+', '-']:
        return strand
    else:
        return {'>': '+', '<': '-'}[strand]


def convertStrandToVista(strand):
    if strand in ['>', '<']:
        return strand
    else:
        return {'>': '+', '<': '-'}[strand]


class Gene(AbstractFeature):
    """Vista gene feature"""
    
    def __init__(self, name='', start=0, end=0, strand=''):
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.strand = convertStrandToPM(strand)
        self.exons = []
        self.flipped = False
    
    def add(self, exon):
        exon.name = self.name
        exon.strand = self.strand
        self.exons.append(exon)
    
    def flip(self, L):
        self.flipped = not self.flipped
        self.start,self.end = L+1-self.end,L+1-self.start
        for exon in self.exons:
            exon.start,exon.end = L+1-exon.end,L+1-exon.start
        self.exons.reverse()
        if self.strand=='+':
            self.strand = '-'
        else:
            self.strand = '+'
    
    def __iter__(self):
        return iter(self.exons)
    
    def __repr__(self):
        output = ['%s %i %i %s' % (self.name, self.start, self.end, self.strand)]
        for exon in self.exons:
            output.append('  %s' % str(exon))
        output.append('\n')
        return '\n'.join(output)
    
    def __getitem__(self, i):
        return self.exons[i]
    
    def __len__(self):
        return len(self.exons)
    
    @staticmethod
    def fromTokens(tokens):
        return Gene(tokens[3], tokens[1], tokens[2], tokens[0])


class Exon:
    def __init__(self, start, end, kind='exon'):
        self.start = int(start)
        self.end = int(end)
        self.kind = kind
    
    def __repr__(self):
        return '%i %i %s' % (self.start, self.end, self.kind)
    
    @staticmethod
    def fromTokens(tokens):
        return Exon(tokens[0], tokens[1], tokens[2])


def load_iter(iFileHandle):
    iFile = smartopen(iFileHandle)
    first = True
    for line in iFile:
        tokens = line.strip().split()
        if not line or line[0]=='#':
            continue
        elif line[0] in ['<', '>']:
            if not first: 
                yield gene
            else:
                first = False
            gene = Gene.fromTokens(tokens)
        else:
            gene.add(Exon.fromTokens(tokens))
    yield gene


def load(iFileHandle):
    genes = []
    for gene in load_iter(iFileHandle):
        genes.append(gene)
    return genes


def load2(iFileHandle):
    iFile = smartopen(iFileHandle)
    genes = []
    for line in iFile:
        tokens = line.strip().split()
        if not line or line[0]=='#':
            continue
        elif line[0] in ['<', '>']:
            genes.append(Gene.fromTokens(tokens))
        else:
            genes[-1].add(Exon.fromTokens(tokens))

