"""
Useful functions and classes
"""

import os, sys, re
from itertools import izip
import warnings


def smartopen(fileHandle, mode='r'):
    """Return an open file object, if passed a filename;
    otherwise just return the object.
    
    @param fileHandle: Input or output file or filename
    @param mode: IO mode ('r', 'w', ...)
    @return: the open file object
    """
    if type(fileHandle)==str:
        ioFile = open(fileHandle, mode)
        return ioFile
    else:
        return fileHandle


def extractRootName(filename):
    """Extract the root name from the filename (path)."""
    return os.path.splitext(os.path.split(filename)[1])[0]


def progressMessage(format, i, n=None, w=80, ostream=sys.stderr):
    if n is None:
        s = ' %i' % i
    else:
        f = 100*float(i)/n
        s = ' %i/%i (%0.1f%%)' % (i,n,f)
    message = format % s
    ostream.write("\b"*w + message)
    ostream.flush()


class ClassFromDict:
    """Simple class that class-ifies a dict"""
    
    def __init__(self, d):
        self.__dict__ = d
    
    def __repr__(self):
        return str(self.__dict__)


def splitBlockCoords(block):
    regex = re.compile('(?P<chrom>\w+)\.(?P<start>\d+)-(?P<end>\d+)')
    rs = regex.search(block)
    d = rs.groupdict()
    return d['chrom'], int(d['start']), int(d['end'])


def unique(a):
    """Return the unique elements of an array."""
    d = {}
    for i,item in enumerate(a):
        try:
            d[item].append(i)
        except KeyError:
            d[item] = [i]
    d = d.items()
    d.sort(key=lambda x: x[1][0])
    return [x[0] for x in d]


def nonunique(x):
    """Return the non-unique elements of an array. The array elements
    must be hashable."""
    d = {}
    for item in x:
        try:
            d[item] += 1
        except KeyError:
            d[item] = 1
    return [item for item,count in d.iteritems() if count>1]


def unique_array_of_dicts(array):
    """Return the unique elements of an array of dictionaries."""
    new = [tuple(d.items()) for d in array]
    return [dict(list(x)) for x in unique(new)]


def unzip(array):
    """Unzips arrays [(a1,b1), (a2,b2), ...]."""
    a = []
    b = []
    for x,y in array:
        a.append(x)
        b.append(y)
    return a,b


def argmax(*args):
    if type(args[0]) in [list, tuple]:
        array = args[0]
    else:
        array = args
    return max(izip(array, xrange(len(array))))


def isort(array):
    """Return sorted array & sorting index."""
    decorated = [(x,i) for i,x in enumerate(array)]
    decorated.sort()
    return unzip(decorated)


def order(array, index):
    """Order array by index provided."""
    return [array[i] for i in index]


def DSU(array, accessor):
    """Decorate Sort Undecorate -- a Schwartzian transform.

    @param array: List to sort.
    @type array: List
    @param accessor: A function of each list element, 
    returning the value by which array will be sorted.    
    e.g. accessor = lambda x: int(x[3]).
    e.g. accessor = lambda x: x.sStart.
    @return: The sorted list.
    """
    decorated = [(accessor(x), x) for x in array]
    decorated.sort()
    return [y[1] for y in decorated]


def loadChrSizes(arg, converter=None):
    """Load chromosome, scaffold or contig sizes.
    
    @param arg: filename or config object
    @param converter: chromosome name converter function (optional)
    """
    
    filename = ''
    if type(arg)==str: # filename
        filename = arg
    else: # config object
        filename = arg.chrSizes
    
    iFile = open(filename)
    chrSizes = {}
    for line in iFile:
        line = line.replace('\t', ' ')
        chrom,size = line.strip().split()
        if converter:
            chrom = converter(chrom)
        chrSizes[chrom] = int(size)
    return chrSizes


def warnDeprecated(message):
    """Issue a deprecation warning. "Deprecated: $message"
    
    @param message: Typically "sequence.reverse_complement (use reverseComplement)."
    """
    warnings.warn(message, DeprecationWarning)    


def raiseDeprecated(message):
    """Raise a deprecation exception.
    @param message: Typically "sequence.reverse_complement (use reverseComplement)."
    """
    raise DeprecationWarning(message)
