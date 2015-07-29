"""
align module
"""

import os, sys, copy
import math, numpy
from fasta import FastaFile
import sequence
from useful import extractRootName, smartopen


# Utility function for loaders
def firstWord(header):
    return header.split()[0]


class Alignment(object):
    def __init__(self, order=None, seq=None, gapChar='-'):
        if not order and not seq:
            self.seqDict = {}
            self.order = []
        else:
            self.seqDict = dict(zip(order, seq))
            self.order = order
        self.gapChar = gapChar
    
    def __len__(self):
        return len(self.seqDict.values()[0])
    
    def numberOfSeqs(self):
        return len(self.seqDict)
    
    def add(self, name, seq):
        if not name in self.seqDict:
            self.seqDict[name] = seq
            self.order.append(name)
        else:
            raise Exception('Sequence already added.')
    
    def append(self, name, seq, cleanup=False):
        if cleanup: 
            seq = ''.join(seq.split())
        
        if not name in self.seqDict:
            self.seqDict[name] = seq
            self.order.append(name)
        else:
            self.seqDict[name] += seq
    
    def remove(self, name):
        del self.seqDict[name]
        self.order.remove(name)
    
    def getName(self, i):
        return self.order[i]
    
    def getSeq(self, i):
        return self.seqDict[self.order[i]]
    
    def getItems(self, i):
        name = self.order[i]
        return name, self.seqDict[name]
    
    def __getitem__(self, name):
        return self.seqDict[name]
    
    def __setitem__(self, name, seq):
        if not name in self.order:
            self.order.append(name)
        self.seqDict[name] = seq
    
    def rename(self, name, newName):
        i = self.order.index(name)
        seq = self.seqDict[name]
        self.order[i] = newName
        self.seqDict[newName] = seq
        del self.seqDict[name]
    
    def __repr__(self):
        output = []
        for name in self.order:
            output.append('%-19s %s' % (name, self.seqDict[name]))
        return '\n'.join(output)
    
    def __iter__(self):
        self.index = 0
        return self
    
    def next(self):
        if self.index==len(self.order):
            raise StopIteration
        self.index += 1
        return self.order[self.index-1]
    
    def column(self, i):
        return [self.seqDict[name][i] for name in self.order]
    
    def column_iter(self):
        for i in xrange(self.__len__()):
            yield self.column(i)
    
    def addColumn(self, col):
        for name,x in zip(self.order, col):
            self.seqDict[name] = ''.join([self.seqDict[name], x])
    
    def deleteColumn(self, i):
        seqDict2 = {}
        for name in self.order:
            seq = list(self.seqDict[name])
            del seq[i]
            seqDict2[name] = ''.join(seq)
        self.seqDict = seqDict2
    
    def columnFreq(self, i, **kw):
        return sequence.calcFrequency(self.column(i), **kw)
    
    def columnFreq_iter(self, **kw):
        for i in xrange(self.__len__()):
            yield self.columnFreq(i, **kw)
    
    def getSubAlignment(self, start, end):
        seqs = []
        for spp in self.order:
            seqs.append(self.seqDict[spp][start:end])
        return Alignment(self.order, seqs)
    
    @staticmethod
    def load(iFileHandle, format='fasta', **kw):
        if format=='fasta':
            align = Alignment.loadFasta(iFileHandle, **kw)
        elif format in ['clustal', 'clustalw']:
            align = Alignment.loadClustal(iFileHandle, **kw)
        elif format=='phylip':
            align = Alignment.loadPhylip(iFileHandle, **kw)
        elif format=='stockholm':
            align = Alignment.loadStockholm(iFileHandle, **kw)
        else:
            raise Exception('Unsupported file format.')
        return align
    
    def save(self, filename, format='fasta', **kw):
        if format=='fasta':
            self.saveFasta(filename, **kw)
        elif format in ['clustal', 'clustalw']:
            self.saveClustal(filename, **kw)
        elif format=='phylip':
            self.savePhylip(filename, **kw)
        elif format=='molphy':
            self.saveMolphy(filename, **kw)
        elif format=='nexus':
            self.saveNexus(filename, **kw)
        else:
            raise Exception('Unsupported file format.')
    
    def saveFasta(self, oFileHandle, **kw):
        writer = FastaFile(oFileHandle, 'w')
        for name in self.order:
            writer.write(name, self.seqDict[name])
        writer.close()
    
    def saveClustal(self, oFileHandle, nameWidth=None, width=60, interleaved=True, **kw):
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, 'CLUSTAL W (1.83) multiple sequence alignment\n\n'
        
        if not nameWidth:
            nameWidth = max([len(name) for name in self.order])
        format = '%%-%is %%s' % nameWidth
        
        L = self.__len__()
        if interleaved:
            for i in xrange(0,L,width):
                for name in self.order:
                    print >> oFile, format % (name[0:nameWidth], self.seqDict[name][i:i+width])
                print >> oFile, format % (' '*nameWidth, ' '*width)
                print >> oFile
        else:
            for name in self.order:
                print >> oFile, format % (name[0:nameWidth], self.seqDict[name])
            print >> oFile, format % (' '*nameWidth, ' '*len(self.seqDict[name]))
            print >> oFile
        oFile.close()
    
    def savePhylip(self, oFileHandle, width=10, **kw):
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, '%7i%7i' % (self.numberOfSeqs(), self.__len__())
        
        L = self.__len__()
        format = '%%-%is %%s' % width
        for i in xrange(0,L,50):
            for name in self.order:
                if i==0:
                    label = name[0:width]
                else:
                    label = ''
                
                print >> oFile, ('%%-%is' % width) % label,
                for j in xrange(0,50,10):
                    print >> oFile, self.seqDict[name][i+j:i+j+10],
                print >> oFile
            print >> oFile
        oFile.close()
    
    def saveMolphy(self, oFileHandle, width=60):
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, '%i %i' % (self.numberOfSeqs(), self.__len__())
        
        for name in self.order:
            print >> oFile, name
            for i in xrange(0, self.__len__(), width):
                print >> oFile, self.seqDict[name][i:i+width]
        oFile.close()
    
    def saveNexus0(self, oFileHandle, datatype='protein', gap='-', interleave=False, width=60):
        assert datatype.lower() in ['dna', 'rna', 'protein', 'standard', 'restriction']
        yesno = {True: 'yes', False: 'no'}
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, '#NEXUS'
        print >> oFile, 'BEGIN data;'
        print >> oFile, 'DIMENSIONS ntax=%i nchar=%i;' % (self.numberOfSeqs(), self.__len__())
        print >> oFile, 'FORMAT datatype=%s interleave=%s gap=%s;' \
            % (datatype, yesno[interleave], gap)
        print >> oFile
        
        print >> oFile, 'MATRIX'
        format = '%%-%is %%s' % max([len(name) for name in self.seqDict])
        if interleave:
            for i in xrange(0,self.__len__(),width):
                for name in self.order:
                    print >> oFile, format % (name, self.seqDict[name][i:i+width])
                print >> oFile
        else:
            for name in self.order:
                print >> oFile, format % (name, self.seqDict[name])
            print >> oFile
        print >> oFile, ';'
        
        print >> oFile, 'END;'
        oFile.close()
    
    def saveNexus(self, oFileHandle, datatype='protein', gap='-', interleave=False, width=60):
        assert datatype.lower() in ['dna', 'rna', 'protein', 'standard', 'restriction']
        yesno = {True: 'yes', False: 'no'}
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, '#nexus'
        print >> oFile
        print >> oFile, 'BEGIN Taxa;'
        print >> oFile, 'DIMENSIONS ntax=%i;' % self.numberOfSeqs()
        print >> oFile, 'TAXLABELS'
        for i,name in enumerate(self.order):
            print >> oFile, "[%i] '%s'" % (i+1,name)
        print >> oFile, ';'
        print >> oFile, 'END; [Taxa]'
        print >> oFile
        
        print >> oFile, 'BEGIN Characters;'
        print >> oFile, 'DIMENSIONS nchar=%i;' % self.__len__()
        print >> oFile, 'FORMAT'
        print >> oFile, '        datatype=%s' % datatype
        print >> oFile, '        missing=?'
        print >> oFile, '        gap=%s' % gap
        print >> oFile, '        symbols="a r n d c q e g h i l k m f p s t w y v z"'
        print >> oFile, '        labels=left'
        print >> oFile, '        transpose=no'
        print >> oFile, '        interleave=%s' % yesno[interleave]
        print >> oFile, ';'
        
        print >> oFile, 'MATRIX'
        format = "%%-%is %%s" % max([len(name) for name in self.seqDict])
        if interleave:
            for i in xrange(0,self.__len__(),width):
                for name in self.order:
                    print >> oFile, format % (name, self.seqDict[name][i:i+width])
                print >> oFile
        else:
            for name in self.order:
                print >> oFile, format % (name, self.seqDict[name])
            print >> oFile
        print >> oFile, ';'
        
        print >> oFile, 'END;'
        oFile.close()
    
    @staticmethod
    def loadFasta(iFileHandle, renamer=firstWord):
        """
        Load fasta alignment data.
        """
        aln = Alignment()
        for h,s in FastaFile(iFileHandle):
            try:
                aln.add(renamer(h), s)
            except Exception, e:
                print e
        return aln
    
    @staticmethod
    def loadPecan(iFileHandle, renamer=firstWord):
        """
        Load pecan alignment data.
        """
        aln = Alignment()
        for h,s in FastaFile(iFileHandle):
            try:
                aln.add(renamer(h), s)
            except Exception, e:
                print e
        
        refSeq = list(aln[aln.order[0]])
        for spp in aln.order[1:]:
            seq = list(aln[spp])
            for i in xrange(len(seq)):
                if seq[i]=='.':
                    seq[i] = refSeq[i]
            aln[spp] = ''.join(seq)
        return aln
    
    @staticmethod
    def loadClustal(iFileHandle, headerCheck=True):
        """
        Load clustal alignment data.
    
        @param iFileHandle: Input filename or file.
        @param headerCheck: Test that CLUSTAL appears on first line (default True)
        """
        iFile = smartopen(iFileHandle)
        clustalHeader = iFile.readline().strip()
        if headerCheck and not clustalHeader.split()[0] in ['CLUSTAL', 'MUSCLE']:
            raise 'Not a CLUSTAL file'
    
        alignment = Alignment()
        for line in iFile:
            if line[0]!=' ':
                tokens = line.strip().split()
                if len(tokens)==2:
                    alignment.append(tokens[0], tokens[1])
        iFile.close()
        return alignment
    
    @staticmethod
    def loadPhylip(iFileHandle, multipleDatasets=False):
        """
        Load phylip alignment data.
        
        @param iFileHandle: Input filename or file.
        @param multipleDatasets: Do not close file (default False)
        """
        iFile = smartopen(iFileHandle)
        nSeq,L = [int(x) for x in iFile.readline().strip().split()]
        
        alignment = Alignment()
        for i in xrange(nSeq):
            line = iFile.readline()
            tokens = line.strip().split()
            seq = ''.join(tokens[1:])
            alignment.append(tokens[0], seq)
        skip = iFile.readline()
        
        width = len(alignment)
        nBlocks = int(math.ceil(float(L)/width))
        for j in xrange(nBlocks-1):
            for i in xrange(nSeq):
                line = iFile.readline().strip()
                name = alignment.getName(i)
                alignment.append(name, line, cleanup=True)
            
            try:
                skip = iFile.readline()
            except:
                pass
        
        if not multipleDatasets: iFile.close()
        
        return alignment
    
    @staticmethod
    def loadStockholm(iFileHandle, **kw):
        """Load Stockholm alignment data.
        
        @param iFileHandle: Input filename or file.
        @returns: a dictionary of sequences {name1: seq1, name2: seq2, ...}
        """
        iFile = smartopen(iFileHandle)
        alignment = Alignment()
        alignment.headers = {}
        
        # Skip hmmalign header
        start = False
        for line in iFile:
            line = line.strip()
            if line=='# STOCKHOLM 1.0':
                break
        
        # Parse sto header info
        for line in iFile:
            line = line.strip()
            if not line: # Blank lines
                continue
            elif line[0:4]=='#=GS': # Fasta headers
                header = line.strip()[1:]
                tokens = header.split()
                alignment.headers[tokens[0]] = header
                continue
            elif line[0]=='#': # Other boring comment lines
                continue
            elif line=='//': # End of file
                break
            
            # The real stuff
            name,seq = line.split()
            alignment.append(name, seq)
        return alignment
    
    @staticmethod
    def fromColumns(species, columns):
        aln = Alignment()
        aln.order = species
        for column in columns:
            for spp,aa in zip(species, column):
                try:
                    aln[spp].append(aa)
                except KeyError:
                    aln[spp] = [aa]
        for spp in aln.seqDict:
            aln.seqDict[spp] = ''.join(aln.seqDict[spp])
        return aln


class Transform:
    """Transform from gapped to ungapped coordinates."""
    def __init__(self, gappedSeq, gapChar='-'):
        """
        @param gappedSeq: Gapped sequence alignment
        """
        self.gappedSeq = gappedSeq
        self.ungappedSeq = ''.join(gappedSeq.split(gapChar))
        self.gapChar = gapChar
        self.g2u = []
        self.u2g = []
        
        u = 0
        g = 0
        for symbol in list(self.gappedSeq):
            self.g2u.append(u)
            if symbol!=gapChar:
                self.u2g.append(g)
                u += 1
            g += 1
    
    def getUngappedLength(self):
        return len(self.ungappedSeq)
    
    def getGappedLength(self):
        return len(self.gappedSeq)
    
    def toGapped(self, i, flip=False):
        """
        Convert ungapped coordinate to gapped.
        @param i: Ungapped coordinate (indexed from 0)
        """
        try:
            if flip:
                L = len(self.gappedSeq)+1
                return L-self.u2g[i]
            else:
                return self.u2g[i]
        except Exception, e:
            print i, len(self.u2g)
            print e
            import traceback
            traceback.print_exc()
            sys.exit()
    
    def toUngapped(self, i, flip=False):
        """
        Convert gapped coordinate to ungapped.
        @param i: Gapped coordinate (indexed from 0)
        """
        if flip:
            L = len(self.ungappedSeq)+1
            return L-self.g2u[i]
        else:
            return self.g2u[i]


def makeReferencedAlignment(seq1, seq2, gapChar='-'):
    oseq1 = []
    oseq2 = []
    for a,b in zip(list(seq1), list(seq2)):
        if a!=gapChar:
            oseq1.append(a)
            oseq2.append(b)
    oseq1 = ''.join(oseq1)
    oseq2 = ''.join(oseq2)
    return oseq1, oseq2


def getUngappedBlocks(gappedSeq, gapChar='-'):
    ungappedSeqBlocks = gappedSeq.split(gapChar)
    ungappedBlocks = []
    pos = 0
    for ungappedSeqBlock in ungappedSeqBlocks:
        L = len(ungappedSeqBlock)
        if L==0:
            pos += 1
        else:
            start = pos
            end = pos+L-1
            ungappedBlocks.append([start, end])
            pos += L+1
    return ungappedBlocks


def calcPercentIdInUngappedBlocks(seq1, seq2, gapChar='-', **kw):
    i = 0
    identityCount = [0]
    blocks = [[]]
    for x,y in zip(list(seq1),list(seq2)):
        if 'N' in [x,y]:
            pass
        elif x!=gapChar and y!=gapChar:
            blocks[-1].append(i)
            if x==y: 
                identityCount[-1] += 1
        else:
            identityCount.append(0)
            blocks.append([])
        i += 1
    
    pids = []
    blocks2 = []
    for b,c in zip(blocks,identityCount):
        if len(b)>=2:
            start = b[0]
            end = b[-1]
            pid = 100.0*c/(end-start+1)
            blocks2.append([start, end, pid])
    return blocks2


def buildIdentityVector(seq1, seq2, gapChar='-'):
    """
    Calculate the percent identity for a pair of sequences from a pairwise or
    multiple alignment.
    """
    identity = numpy.zeros(len(seq1), dtype=float)
    i = 0
    for x,y in zip(list(seq1),list(seq2)):
        if not 'N' in [x,y]:
            pass
        elif x!=gapChar and x==y:
            identity[i] = 100.0
        i += 1
    return identity


def buildIdentityVector2(seq1, seq2, gapChar='-'):
    """
    Calculate the percent identity for a pair of sequences from a pairwise or
    multiple alignment.
    """
    identity = numpy.zeros(len(seq1), dtype=float)
    i = 0
    for x,y in zip(list(seq1),list(seq2)):
        if gapChar in [x,y]:
            identity[i] = -1
        elif x==y:
            identity[i] = 100.0
        i += 1
    return identity


def downSampleToBlocks(data, n=500):
    L = len(data)
    w = L/n
    scale = 1.0/w
    vector = numpy.array(data, dtype=float)
    blocks = []
    for i in xrange(n):
        start = i*w
        if i==n-1:
            end = L-1
        else:
            end = (i+1)*w
        pid = scale*numpy.sum(vector[start:end+1])
        blocks.append([start, end, pid])
    return blocks


def downSampleToBlocks2(data, n=500):
    L = len(data)
    w = L/n
    vector = numpy.array(data, dtype=float)
    blocks = []
    for i in xrange(n):
        start = i*w
        if i==n-1:
            end = L-1
        else:
            end = (i+1)*w
        mask = numpy.where(vector[start:end+1]!=-1,1,0)
        pid = 100.0*float(numpy.sum(mask))/len(mask)
        blocks.append([start, end, pid])
    return blocks


def calcPercentIdInWindows(seq1, seq2, n=500, gapChar='-'):
    identity = buildIdentityVector2(seq1, seq2, gapChar=gapChar)
    return downSampleToBlocks2(identity, n=n)
