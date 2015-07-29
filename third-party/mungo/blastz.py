"""
blastz module

Example blastz stanza::

    a {
      s 10978
      b 529 48927125
      e 726 48927322
      l 529 48927125 726 48927322 76
    }

"""

import sys, copy
from mungoCore import *
from useful import smartopen


strandDict = {'0': '+', '1': '-'}


class HSP:
    def __init__(self, tokens):
        self.interval1 = [int(tokens[0]), int(tokens[2])]
        self.interval2 = [int(tokens[1]), int(tokens[3])]
        self.pcIdent = int(tokens[4])
    
    def __repr__(self):
        return '%i %i %i %i %i' % (self.interval1[0], self.interval1[1], 
            self.interval2[0], self.interval2[0], self.pcIdent)


class Chain:
    def __init__(self):
        self.filenames = []
        self.headers = []
        self.lengths = []
        self.strands = []
        
        self.score = 0
        self.interval1 = []
        self.interval2 = []
        self.hsps = []
    
    def __repr__(self):
        return '%s%s-->%s%s Score=%s' \
            % tuple([str(x) for x in [self.interval1, self.strands[0], self.interval2, self.strands[1], self.score]])
    
    @staticmethod
    def fromBlock(block):
        chain = Chain()
        chain.filenames = block.filenames
        chain.headers = block.headers
        chain.lengths = block.lengths
        chain.strands = block.strands
        
        chain.score = 0
        chain.interval1 = [0,0]
        chain.interval2 = [0,0]
        chain.hsps = []
        return chain


class Block:
    def __init__(self):
        self.filenames = []
        self.headers = []
        self.lengths = []
        self.strands = []


class Alignment:
    def __init__(self):
        self.matrix = []
        self.chains = []
    
    def __repr__(self):
        output = []
        for chain in self.chains:
            output.append(str(chain))
        return '\n'.join(output)
    
    @staticmethod
    def load(iFileHandle):
        iFile = smartopen(iFileHandle)
        aln = Alignment()
        state = None
        for line in iFile:
            line = line.rstrip()
            if not line:
                continue
            elif line=='#:lav': # Section break
                state = 'BLOCK'
                block = Block()
                continue
            elif line=='#:eof': # End of file
                state = 'EOF'
                break
            elif line[0]=='d': # Substitution matrix stanza
                state = 'MATRIX'
                continue
            elif line[0]=='s': # Sequence files stanza
                state = 'FILES'
                continue
            elif line[0]=='h': # Fasta headers stanza
                state = 'HEADERS'
                continue
            elif line[0]=='a': # Alignment stanza
                state = 'ALIGN'
                chain = Chain.fromBlock(copy.copy(block))
                continue
            elif line[0] in ['x', 'm']:
                state = 'BORING'
                continue
            elif state=='MATRIX' and line[0]=='}':
                state = 'MATRIX_END'
                aln.matrix = '\n'.join(aln.matrix)
                continue
            elif state=='ALIGN' and line[0]=='}':
                state = 'ALIGN_END'
                aln.chains.append(chain)
                chain = None
                continue
            elif line[0]=='}': # End of state
                continue
            
            tokens = line.lstrip().split()
            if state==None:
                print line
                raise Exception('Wrong')
            elif state=='MATRIX':
                aln.matrix.append(line)
            elif state=='FILES':
                block.filenames.append(tokens[0])
                block.lengths.append(int(tokens[2]))
                block.strands.append(strandDict[tokens[3]])
                # Next line in stanza
                tokens = iFile.next().strip().split()
                block.filenames.append(tokens[0])
                block.lengths.append(int(tokens[2]))
                block.strands.append(strandDict[tokens[3]])
            elif state=='HEADERS':
                block.headers.append(tokens[0])
                # Next line in stanza
                tokens = iFile.next().strip().split()
                block.headers.append(tokens[0])
            elif state=='ALIGN':
                if tokens[0]=='s':
                    chain.score = int(tokens[1])
                elif tokens[0]=='b':
                    chain.interval1[0] = int(tokens[1])
                    chain.interval2[0] = int(tokens[2])
                elif tokens[0]=='e':
                    chain.interval1[1] = int(tokens[1])
                    chain.interval2[1] = int(tokens[2])
                elif tokens[0]=='l':
                    chain.hsps.append(HSP(tokens[1:]))
        return aln


if __name__=='__main__':
    aln = Alignment.load('/Users/papenfuss/tammar/margaret/blastz/chrX.txt')
    i = 0
    for chain in aln.chains:
        i += 1
        print i,chain
