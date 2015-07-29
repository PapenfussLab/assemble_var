"""
cap3 module
"""

import os.path
import re
from align import Alignment
from mungoCore import *
from useful import smartopen
import numpy


class State:
    begin = -1
    beginSection2 = -2
    section1 = 1
    section2 = 2
    section1or2 = 3


class Tokens:
    endSection1 = "DETAILED DISPLAY OF CONTIGS"
    ruler = '                          .    :    .    :    .    :    .    :    .    :    .    :'


def CapFile(iFileHandle, **kw):
    """Factory function for CAP3 Reader classes
    
    @param iFileHandle: CAP3 file name or object
    @keyword force: Force indexing CAP file (Default: False)
    """
    return CapReaderIndexed(iFileHandle, **kw)


class CapReaderIndexed(AbstractDataReader):
    """Class for accessing CAP3 files using an index"""
    
    def __init__(self, iFileHandle, force=False, section=None, **kw):
        """Constructor
        
        @param iFileHandle: CAP3 file name or object
        """
        self.iFile = smartopen(iFileHandle)
        self.iFilename = self.iFile.name
        
        self.indexFile = CapIndexFile(self.iFilename)
        self.indexFile.build(force=force)
        self.stopAtMiddle = False
        if section==1:
            self.stopAtMiddle = True
        elif section==2:
            self.seek(0,2)
        
        self._iter = None
        self._initIter = True
        self._section = section
        
    def reset(self):
        self.iFile.seek(0)
        self._section = 1
        self._iter = self._generator()
    
    def __iter__(self):
        if self._initIter:
            self._iter = self._generator()
            self._initIter = False
        return self
    
    def next(self):
        for header,contents in self._iter:
            contents = parseContents(contents, self._section)
            return header,contents
        raise StopIteration()
    
    def readOne(self):
        if self._initIter:
            self._iter = self._generator()
            self._initIter = False
        return self.next()
    
    def seek(self, i, section=2):
        """Seek position in CAP file
        
        @param i: ith contig (Indexed from 0)
        @keyword section: Section 1 or 2 (Default: 2)
        """
        assert section in (1,2)
        self._section = section
        pos = self.indexFile.get(i)
        if not pos is None:
            self.iFile.seek(pos[section-1])
            self._iter = self._generator()
        else:
            raise EOFError("cap3.CapReaderIndexed.seek")
    
    def get(self, i, section=2):
        """Return slice of CAP file.
        
        @param i: Index or slice
        """
        self._section = section
        if type(i) in (int,long):
            self.seek(i, section=section)
            return self.readOne()
        elif type(i)==slice:
            assert i.step in (None,1)
            self.seek(i.start, section=section)
            j = i.start
            results = []
            for header,contents in self: # self._iter:
                results.append([header,contents])
                j += 1
                if j==i.stop: break
        return results
    
    def search(self, accession, section=2):
        """Search for accession - private mapping interface
        
        @param accession: Accession key
        @keyword section: Section 1 or 2 (Default: 2)
        """
        assert section in (1,2)
        self._section = section
        pos = self.indexFile.search(accession)
        self.iFile.seek(pos[section-1])
        self._iter = self._generator()
        return self.readOne()
    
    def __getitem__(self, key):
        """Define a container-like interface to Section 2 entries."""
        self._section = 2
        if type(key) in(int, slice):
            return self.get(key, section=2)
        elif type(key)==str:
            return self.search(key, section=2)
        else:
            raise Exception('cap3.CapReaderIndexed.__getitem__')
    
    def __len__(self):
        return len(self.indexFile)
    
    def _generator(self):
        header = ''
        contents = []
        state = State.begin
        for line in self.iFile:
            line = line.rstrip()
            if not line:
                # Blank
                pass
            elif state==State.begin:
                if line[0]!='*':
                    # Junk at top
                    pass
                elif line[0]=='*':
                    # First contig
                    state = State.section1or2
                    header = ''.join(line.replace('*', '').strip().split())
            elif state==State.section1or2:
                if line[0]=='*':
                    # New contig
                    yield header, contents
                    header = ''.join(line.replace('*', '').strip().split())
                    contents = []
                elif line==Tokens.endSection1:
                    if self.stopAtMiddle:
                        yield header, contents
                        raise StopIteration()
                    self._section = 2
                else:
                    # Contents line
                    contents.append(line)
        yield header,contents


class CapIndexFile:
    """CAP3 output index file class"""
    
    def __init__(self, iFilename):
        """Constructor
        
        @param iFilename: CAP3 filename
        """
        self.iFilename = iFilename
        self.idxFilename = '%s.index' % iFilename
        self.isCached = False
        self.index = None
    
    def isIndexed(self):
        """Test if index file exists."""
        return os.path.exists(self.idxFilename)
    
    def build(self, force=False):
        """Build index file.
        
        @keyword force: Overwrite existing index file (default=False) 
        """
        if not self.isIndexed() or force:
            idxFile = open(self.idxFilename, 'w')
            order,positions = self._calcPositions()
            for contigName in order:
                pos1,pos2 = positions[contigName]
                idxFile.write("%s\t%i\t%i\n" % (contigName, pos1, pos2))
            idxFile.write('Section1\t%i\t%i\n' % (positions[order[0]][0], positions[order[-1]][0]))
            idxFile.write('Section2\t%i\t%i\n' % (positions[order[0]][1], positions[order[-1]][1]))
            idxFile.close()
    
    def search(self, contigName):
        """Mapping interface.
        
        @param contigName: Contig name
        """
        idxFile = open(self.idxFilename)
        for line in idxFile:
            tokens = line.strip().split('\t')
            name,pos1,pos2 = tokens[0],int(tokens[1]),int(tokens[2])
            if name==contigName:
                return pos1,pos2
    
    def get(self, index):
        """Container interface.
        
        @param index: Index of CAP3 entry (indexed from 0)
        """
        idxFile = open(self.idxFilename)
        j = 0
        for line in idxFile:
            tokens = line.strip().split('\t')
            name,pos1,pos2 = tokens[0],int(tokens[1]),int(tokens[2])
            if j==index:
                return pos1,pos2
            j += 1
    
    def cache(self):
        """Cache index in dict."""
        idxFile = open(self.idxFilename)
        self.index = {}
        for line in idxFile:
            tokens = line.strip().split('\t')
            contigName,pos1,pos2 = tokens[0],int(tokens[1]),int(tokens[2])
            self.index[accession] = (pos1,pos2)
        self.isCached = True
        self.search = self._searchCached
        self.get = self._getCached
    
    def _searchCached(self, accession):
        """Mapping interface when cached."""
        return self.index[accession]
    
    def _getCached(self, index):
        idxList = self.index.items()
        idxList.sort(key=lambda x: x[1])
        return idxList[index][1:]
    
    def __len__(self):
        """Index file length."""
        idxFile = open(self.idxFilename)
        count = 0
        for line in idxFile:
            if line[0]=='>':
                count += 1
        return count
    
    def __getitem__(self, key):
        """Mapping/container interface"""
        if type(key) in [int, long]:
            return self.get(key)
        else:
            return self.search(key)
    
    def _calcPositions(self):
        """Returns index data for a CAP3 file.
        
        @returns: Iterator yielding (name, pos)
        """
        iFile = open(self.iFilename)
        order = []
        positions = {}
        pos = 0
        for line in iFile:
            if line[0]=='*':
                line1 = line.replace('*', '')
                tokens =line1.strip().split()
                contigName = tokens[0]+tokens[1]
                try:
                    positions[contigName].append(pos)
                except KeyError:
                    order.append(contigName)
                    positions[contigName] = [pos]
            pos += len(line)
        iFile.close()
        return order, positions


def CapIndexFactory(iFilename, force=False):
    index = CapIndexFile(iFilename)
    index.build(force=force)
    return index


def parseContents(contents, section):
    if section==1:
        readNames = []
        contained = {}
        for line in contents:
            if line[0]!=' ':
                readNames.append(line.strip())
            else:
                tokens = line.strip().split()
                readNames.append(tokens[0])
                contained[tokens[0]] = tokens[-1]
        return readNames, contained
    elif section==2:
        aln = Alignment()
        iSeqStart = -1
        L = 0
        for line in contents:
            if line==Tokens.ruler:
                # Ruler line
                iSeqStart = line.find('.')-4 # = 22 ???
                if iSeqStart!=22: print "!!! iSeqStart:", iSeqStart
            elif line[0]==' ':
                # Separator line
                pass
            else:
                line = line.rstrip()
                name = line[0:iSeqStart].strip()
                seq = line[iSeqStart:]
                if name!='consensus':
                    if not name in aln.seqDict:
                        seq = seq.replace(' ', '.')
                        seq = '.'*L + seq
                    aln.append(name, seq)
                else:
                    aln.append(name, seq)
                    L = len(aln['consensus'])
        
        # Pad short sequences with '.'
        for name in aln:
            aln[name] = aln[name] + '.'*(L-len(aln[name]))
        return aln
    else:
        raise Exception("cap3.parseContents")


def meanDepth(aln):
    depth = []
    for c in aln.column_iter():
        c = ''.join(c).replace('.', '')
        depth.append(len(c)-1)
    return numpy.mean(depth)


def medianDepth(aln):
    depth = []
    for c in aln.column_iter():
        c = ''.join(c).replace('.', '')
        depth.append(len(c)-1)
    return numpy.median(depth)


def avDepth(aln):
    return meanDepth(aln)


if __name__=='__main__':
    iFilename = '/Users/papenfuss/devil/assembly/unclean/tumour_cap3.out'
    f = CapReaderIndexed(iFilename, force=True)
    f.seek(1,2)
    print f.readOne()
    