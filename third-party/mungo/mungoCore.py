"""
mungoCore - abstract base classes and core utilites
"""

import sys
import copy
from useful import smartopen, raiseDeprecated


def attributesToFormat(attributes):
    return '\t'.join(['%%(%s)s' % attr for attr in attributes])


class AbstractFeature(object):
    """Feature base class"""
    
    def __init__(self, *args, **kw):
        self._rawData = None
        self.indefault = kw.pop('indefault', None)
        self.outdefault = kw.pop('outdefault', '.')
        doConvert = kw.pop('doConvert', True)
        
        if not args and not kw:
            self._fromDict({})
            doConvert = False
        elif kw:
            self._fromDict(kw)
        elif type(args[0])==dict:
            self._fromDict(args[0])
        elif type(args[0]) in [list, tuple]:
            self._fromList(args[0])
        elif len(args)==len(attributes):
            self._fromDict(args)
        else:
            try:
                self._fromDict(args[0].__dict__)
            except:
                print attributes
                print args
                print kw 
                raise Exception, "I don't know what else to do."
        
        if doConvert:
            self.convert()
    
    def convert(self):
        if type(self.converters)==dict:
            self.converters = self.converters.items()
        
        for attrib,converter in self.converters:
            try:
                self.__dict__[attrib] = converter(self.__dict__[attrib])
            except TypeError:
                pass
            except ValueError:
                print self.__dict__
                sys.exit(-1)
    
    def __repr__(self):
        return self.format % self.__dict__
    
    def _fromList(self, tokens):
        self._rawData = tokens
        for key,value in zip(self.attributes, tokens):
            if value!=self.outdefault:
                self.__dict__[key] = value
            else:
                self.__dict__[key] = self.indefault
    
    def _fromDict(self, d):
        self._rawData = d
        for attrib in self.attributes:
            self.__dict__[attrib] = d.get(attrib, self.indefault)
    
    def pretty(self, w=12):
        format = '%%-%is\t%%s' % w
        output = []
        for attrib in self.attributes:
            output.append(format % (attrib, str(self.__dict__[attrib])))
        return '\n'.join(output)
    
    def prettyprint(self):
        print self.pretty()
    
    def addAttribute(self, field, value=None, format=None):
        self.attributes.append(field)
        if not value:
            self.__dict__[field] = self.indefault
        else:
            self.__dict__[field] = value
        
        if not format:
            self.format = self.format + '\t%%(%s)s' % field
        else:
            self.format = self.format + format % field


class AbstractDataFile(object):
    def flush(self):
        self.iFile.flush()
    
    def close(self):
        self.iFile.close()


class AbstractDataReader(AbstractDataFile):
    def __init__(self, iFileHandle):
        self.iFile = smartopen(iFileHandle)
        self.iFilename = self.iFile.name
        self._iter = None
    
    def _generator(self):
        pass
    
    def __iter__(self):
        self._iter = self._generator()
        return self
    
    def next(self):
        for x in self._iter:
            return x
        raise StopIteration
    
    def readAll(self):
        """Read all entries.
        
        @returns: A list
        """
        results = []
        for result in self:
            results.append(result)
        return results
