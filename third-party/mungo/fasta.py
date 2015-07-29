"""
fasta module
"""

import os.path
import re
import cPickle
import sqlite3
import tempfile
from mungoCore import *
from useful import smartopen, progressMessage


def FastaFile(iFileHandle, mode='r', indexed=False, **kw):
    """Factory function for Reader and Writer classes.
    
    @param iFileHandle: Fasta file name or object
    @keyword mode: read(r), write(w) or append(a)
    @keyword indexed: Index the fasta file; reading only (default False)
    """
    if mode=='r' and not indexed:
        return FastaReader(iFileHandle, **kw)
    elif mode=='r' and indexed:
        return FastaReaderIndexed(iFileHandle, **kw)
    elif mode in ('w', 'a'):
        return FastaWriter(iFileHandle, mode=mode, **kw)


class FastaWriter(AbstractDataFile):
    """Class for writing fasta files."""
    
    def __init__(self, iFileHandle, mode='w', width=60, blockSize=None, **kw):
        """
        @param iFileHandle: Output file or name
        @keyword mode: File mode - write(w) or append(a)
        """
        assert mode in ('w', 'a')
        self.iFile = smartopen(iFileHandle, mode)
        self.iFilename = self.iFile.name
        self.width = width
        
        if blockSize: 
            self.setBlockSize(blockSize)
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def setBlockSize(self, blockSize):
        """Split sequences into blocks.
        
        @param blockSize: Sequence block size, when used this is typically ~5Mb
        """
        self.blockSize = blockSize
        self.write = self._writeBlock
    
    def write(self, header, seq):
        """Write a fasta entry.
        
        @param header: Fasta header
        @param seq: Sequence
        """
        print >> self.iFile, '>%s' % header
        print >> self.iFile, pretty(seq, width=self.width)
    
    def _writeBlock(self, header, seq):
        """Alternative write function when using blockSizes"""
        tokens = header.split()
        if len(tokens)==1:
            newHeaderTemplate = '>%s:%%i-%%i' % tokens[0]
        else:
            newHeaderTemplate = '>%s:%%i-%%i %s' % (tokens[0], ' '.join(tokens[1:]))
        
        for j in xrange(0,len(seq),self.blockSize):
            s = seq[j:j+self.blockSize]
            size = min(len(s), self.blockSize)
            print >> self.iFile, newHeaderTemplate % (j+1, j+size)
            print >> self.iFile, pretty(s, width=self.width)
    
    def __call__(self, header, seq):
        """Call interface to write."""
        self.write(header, seq)


class FastaReader(AbstractDataReader):
    """Simple class for reading fasta files"""
    
    def __init__(self, iFileHandle, **kw):
        """
        @param iFileHandle: Fasta file name or object
        """
        AbstractDataReader.__init__(self, iFileHandle)
        self._initIter = True
        self._isAsMapping = False
        self._seqDict = None
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def reset(self):
        self.iFile.seek(0)
        self._iter = self._generator()
    
    def __iter__(self):
        if self._initIter:
            self._iter = self._generator()
            self._initIter = False
        return self
    
    def next(self):
        for h,s in self._iter:
            return h,s
        raise StopIteration
    
    def _generator(self):
        header = ''
        seq = []
        for line in self.iFile:
            line = line.strip()
            if line:
                if line[0]=='>':
                    if seq:
                        yield header, ''.join(seq)
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line)
        if header and seq:
            yield header, ''.join(seq)
    
    def lengthGenerator(self):
        header = ''
        L = 0
        for line in self.iFile:
            line = line.strip()
            if line:
                if line[0]=='>':
                    if L:
                        yield header, L
                    header = line[1:]
                    L = 0
                else:
                    L += len(line)
        if header and L:
            yield header, L
    
    def readAll(self):
        """Read all fasta entries.
        
        @returns: A list of 2-ples [(header,sequence), ...]
        """
        self.reset()
        seqs = []
        for h,s in self:
            seqs.append((h,s))
        return seqs
    
    def readOne(self):
        if self._initIter:
            self._iter = self._generator()
            self._initIter = False
        return self.next()


class Interface:
    CONTAINER = 'container'
    MAPPING = 'mapping'

class IndexMethod:
    SQLITE = 'sqlite3'
    PICKLE = 'pickle'
    TEXT = 'text'

class FastaReaderIndexed(FastaReader):
    """Class for accessing fasta files using an index"""
    
    def __init__(self, iFileHandle, clobber=False, 
        interface=Interface.CONTAINER, method=IndexMethod.SQLITE, **kw):
        """
        @param iFileHandle: Fasta file name or object
        """
        self.iFile = smartopen(iFileHandle)
        self.iFilename = self.iFile.name
        if method==IndexMethod.PICKLE:
            self.indexFile = FastaIndexPickleFile(self.iFilename)
        elif method==IndexMethod.TEXT:
            self.indexFile = FastaIndexTextFile(self.iFilename)
        else: # sqlite3 method is the default
            self.indexFile = FastaIndexFile(self.iFilename)
        
        self.indexFile.build(clobber=clobber)
        self.interface = interface
        self._iter = None
        self._initIter = True
        
        if kw: 
            print 'Uncaptured keywords'
            print kw
    
    def seek(self, i):
        """Seek position in fasta file
        
        @param i: ith fasta entry
        """
        pos = self.indexFile.get(i)
        if not pos is None:
            self.iFile.seek(pos)
            self._iter = self._generator()
        else:
            raise EOFError("FastaReaderIndexed.seek")
    
    def get(self, index):
        """Return slice of fasta file.
        
        @param index: Index or slice
        """
        if type(index) in (int,long):
            self.seek(index)
            return self.readOne()
        elif type(index)==slice:
            assert index.step in (None,1)
            self.seek(index.start)
            j = index.start
            seqs = []
            for h,s in self._iter:
                seqs.append((h,s))
                j += 1
                if j==index.stop: break
            return seqs
    
    def search(self, accession):
        """Search for accession - private mapping interface
        
        @param accession: Accession key
        """
        pos = self.indexFile.search(accession)
        self.iFile.seek(pos)
        self._iter = self._generator()
        return self.readOne()
    
    def like(self, accession):
        """Search for accession - private mapping interface
        
        @param accession: Accession key
        """
        pos = self.indexFile.like(accession)
        self.iFile.seek(pos)
        self._iter = self._generator()
        return self.readOne()
    
    def __getitem__(self, index):
        if self.interface==Interface.CONTAINER:
            return self.get(index)
        elif self.interface==Interface.MAPPING:
            return self.search(index)
        else:
            raise Exception('FastaReaderIndexed: No container-type set.')
    
    def asContainer(self):
        """Define a container-like interface to the fasta entries."""
        self.interface = Interface.CONTAINER
    
    def asMapping(self):
        """Define a mapping interface to the fasta entries.
        
        @keyword convert: Anonymous function defining how to convert fasta header into accession
        """
        self.interface = Interface.MAPPING
    
    def __len__(self):
        return len(self.indexFile)


class FastaIndexFile:
    """Fasta index file class"""
    
    def __init__(self, iFilename):
        """
        @param iFilename: Fasta filename
        """
        self.iFilename = iFilename
        self.idxFilename = '%s_idx.sqlite3' % iFilename
        self.connection = sqlite3.connect(self.idxFilename)
        self.connection.isolation_level = None
        self.cursor = self.connection.cursor()
    
    def isIndexed(self):
        """Test if index file exists."""
        try:
            self.cursor.execute('select * from Offsets limit 1;')
            return len(self.cursor.fetchone())!=0
        except:
            return False
    
    def build(self, clobber=False, separator='!'):
        """Build index file.
        
        @keyword clobber: Overwrite existing index file (default=False) 
        """
        if clobber:
            try:
                self.cursor.execute('drop table Offsets;')
            except:
                pass
        
        if not self.isIndexed():
            schema = """
                CREATE TABLE Offsets (
                    id INTEGER,
                    accession TEXT,
                    offset INTEGER
                    );
                    create index idx_offsets_id on Offsets (id);
                    create index idx_offsets_accession on Offsets (accession);
            """
            self.connection.executescript(schema)
            tmpFile = tempfile.NamedTemporaryFile()
            for i,(accession,offset) in enumerate(self._byteOffsetGenerator()):
                print >> tmpFile, "%i%s%s%s%i" % (i,separator,accession,separator,offset)
                if i % 1000==0:
                    progressMessage("Sequences: %s", i+1)
            progressMessage("Sequences: %s\n", i+1)
            tmpFile.flush()
            
            cmd = """sqlite3 -separator '%s' %s '.import "%s" Offsets'""" \
                % (separator, self.idxFilename, tmpFile.name)
            os.system(cmd)
            tmpFile.close()
    
    def search(self, accession):
        """Mapping interface.
        
        @param accession: Accession key
        """
        self.cursor.execute("select offset from Offsets where accession='%s';" % accession)
        result = self.cursor.fetchone()
        if result:
            return result[0]
        else:
            raise KeyError()
    
    def like(self, accession):
        """Mapping interface.
        
        @param accession: Accession key
        """
        self.cursor.execute("select offset from Offsets where accession glob '%s*';" % accession)
        result = self.cursor.fetchone()
        if result:
            return result[0]
        else:
            raise KeyError()
    
    def get(self, index):
        """Container interface.
        
        @param index: Index of fasta entry
        """
        self.cursor.execute("select offset from Offsets where id=%i;" % index)
        result = self.cursor.fetchone()
        if result:
            return result[0]
        else:
            raise IndexError()
    
    def __len__(self):
        """Fasta index length."""
        self.cursor.execute("select count(*) from Offsets;")
        return self.cursor.fetchone()
    
    def __getitem__(self, key):
        """Mapping/container interface"""
        if type(key) in [int, long]:
            return self.get(key)
        else:
            return self.search(key)
    
    def _byteOffsetGenerator(self):
        """Return an iterator to a multi-fasta file.
        
        @returns: Iterator yielding (name, pos)
        """
        iFile = open(self.iFilename)
        pos = 0
        for line in iFile:
            if line[0]=='>':
                name = line[1:].strip().split()[0]
                yield name, pos
            pos += len(line)


class FastaIndexPickleFile:
    """Fasta index file using a pickled list"""
    
    def __init__(self, iFilename):
        """
        @param iFilename: Fasta filename
        """
        self.iFilename = iFilename
        self.idxFilename = '%s_idx.pkl' % iFilename
        self.index = []
    
    def isIndexed(self):
        """Test if index file exists."""
        return os.path.exists(self.idxFilename)
    
    def build(self, clobber=False):
        """Build index file.
        
        @keyword clobber: Overwrite existing index file (default=False) 
        """
        if not self.isIndexed() or clobber:
            idxFile = open(self.idxFilename, 'w')
            for accession,offset in self._byteOffsetGenerator():
                self.index.append((accession,offset))
            cPickle.dump(self.index, open(self.idxFilename, 'w'))
        else:
            self.index = cPickle.load(open(self.idxFilename))
    
    def search(self, accession):
        """Mapping interface.
        
        @param accession: Accession key
        """
        d = dict(self.index)
        return d[accession]
    
    def get(self, index):
        """Container interface.
        
        @param index: Index of fasta entry
        """
        return self.index[index][1]
    
    def __len__(self):
        """Fasta index length."""
        return len(self.index)
    
    def __getitem__(self, key):
        """Mapping/container interface"""
        if type(key) in [int, long]:
            return self.get(key)
        else:
            return self.search(key)
    
    def _byteOffsetGenerator(self):
        """Return an iterator to a multi-fasta file.
        
        @returns: Iterator yielding (name, pos)
        """
        iFile = open(self.iFilename)
        pos = 0
        for line in iFile:
            if line[0]=='>':
                name = line[1:].strip().split()[0]
                yield name, pos
            pos += len(line)


class FastaIndexTextFile:
    """Text file-based fasta index file class"""
    
    def __init__(self, iFilename):
        """
        @param iFilename: Fasta filename
        """
        self.iFilename = iFilename
        self.idxFilename = '%s_idx.txt' % iFilename
    
    def isIndexed(self):
        """Test if index file exists."""
        return os.path.exists(self.idxFilename)
    
    def build(self, clobber=False):
        """Build index file.
        
        @keyword clobber: Overwrite existing index file (default=False) 
        """
        if not self.isIndexed() or clobber:
            idxFile = open(self.idxFilename, 'w')
            for h,offset in self._byteOffsetGenerator():
                idxFile.write("%i\t%s\n" % (offset, h))
            idxFile.close()
    
    def search(self, accession):
        """Mapping interface.
        
        @param accession: Accession key
        """
        idxFile = open(self.idxFilename)
        for line in idxFile:
            offset,name = line.strip().split('\t')
            if name==accession:
                return int(offset)
    
    def get(self, index):
        """Container interface.
        
        @param index: Index of fasta entry
        """
        idxFile = open(self.idxFilename)
        j = 0
        for line in idxFile:
            tokens = line.strip().split('\t')
            pos,name = int(tokens[0]),tokens[1]
            if j==index:
                return pos
            j += 1
        return None
    
    def __len__(self):
        """Fasta index length."""
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
    
    def _byteOffsetGenerator(self):
        """Return an iterator to a multi-fasta file.
        
        @returns: Iterator yielding (name, pos)
        """
        iFile = open(self.iFilename)
        pos = 0
        for line in iFile:
            if line[0]=='>':
                name = line[1:].strip().split()[0]
                yield name, pos
            pos += len(line)


class Sequence:
    """Fasta sequence class"""
    
    def __init__(self, header='', seq=''):
        self.header = header
        self.seq = seq
    
    def __repr__(self):
        out = ">%s\n%s" % (self.header, pretty(self.seq))
        return out


# Utility functions

def FastaIndexFactory(iFilename):
    index = FastaIndexFile(iFilename)
    index.build()
    return index


def pretty(seq, width=60, header=None, joinChar='\n'):
    """
    Return a prettified version of seq. Default returns the string reformatted
    to 60 chars wide, e.g.
    
    pretty(seq, width=10, joinChar=' ') returns the string with a space every 10 chars.
    
    @param seq: Sequence string
    @param width: Sequence width (default 60)
    @param joinChar: Character to join on (default \\n)
    @rtype: string
    @return: a pretty-version of seq.
    """
    output = []
    if header: output.append(">%s" % header)
    for i in xrange(0, len(seq), width):
        output.append(seq[i:i+width])
    return joinChar.join(output)


def extractSpp(header):
    sppPattern = re.compile('\[.*\]')
    spp = sppPattern.findall(h)[-1][1:-1]
    return spp


def fastaCount(iFilename):
    count = 0
    for line in open(iFilename):
        if line[0]=='>':
            count += 1
    return count
