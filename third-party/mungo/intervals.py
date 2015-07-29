"""
interval module - Interval features

To do:
    - Need to add simple classes or functions for dealing with lists of intervals
"""

import copy
from mungoCore import *
import hmmer, blast, fasta
from useful import DSU


class Interval(object):
    attributes = ['chrom', 'start', 'end', 'strand', 'value']
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        self.n = 1
        for attrib in self.attributes:
            self.__dict__[attrib] = None
        
        if len(args)>0:
            for i,arg in enumerate(args):
                self.__dict__[self.attributes[i]] = arg
        
        if len(kw)>0:
            for attrib in kw:
                self.__dict__[attrib] = kw[attrib]
        
        if (not self.strand or self.strand=='+') and self.start>self.end:
            self.start,self.end = self.end,self.start
            self.strand = '-'
    
    @staticmethod
    def fromList(intervalList):
        """
        Static method; creates an Interval object from a list or tuple
        @param intervalList: list or tuple [start, end]
        @returns: Interval object
        """
        f = Interval(intervalList)
        return f
    
    @staticmethod
    def fromObject(obj):
        """
        Static method; creates an interval object from a HSP or Domain object
        
        @param obj: HSP or Domain object
        @returns: Interval object
        """
        f = Interval()
        f.origObject = obj
        f.start = obj.sStart
        f.end = obj.sEnd
        f.value = obj
        if obj.__class__==hmmer.BlockSixFrameDomain:
            f.chrom = obj.accession
            f.strand = obj.strand
        elif obj.__class__==hmmer.SixFrameDomain:
            f.chrom = obj.accession
            f.strand = obj.strand
        elif obj.__class__==blast.HSP:
            f.chrom = obj.subjectId
            f.strand = obj.strand()
        return f
    
    def toList(self):
        return self.chrom, self.start, self.end, self.strand, self.value
    
    def toDict(self):
        return dict(zip(self.attributes, self.toList()))
    
    def __repr__(self):
        return self.format % self.toDict()
    
    def overlaps(self, f):
        status = False
        if self.chrom==f.chrom and self.strand==f.strand:
            if self.start<=f.start<=self.end:
                status = True
            elif self.start<=f.end<=self.end:
                status = True
            elif f.start<=self.start<=f.end:
                status = True
        return status
    
    def strandlessOverlaps(self, f):
        status = False
        if self.chrom==f.chrom:
            if self.start<=f.start<=self.end:
                status = True
            elif self.start<=f.end<=self.end:
                status = True
            elif f.start<=self.start<=f.end:
                status = True
        return status
    
    def __cmp__(self, other):
        return cmp(self.start, other.start) or cmp(self.end, other.end)
    
    def merge(self, f):
        assert self.start<self.end
        assert f.start<f.end
        self.start = min(self.start, f.start)
        self.end = max(self.end, f.end)
        self.n += 1
    
    def getSequence(self, blastDb, padFivePrime=0, padThreePrime=0):
        start = max(1,self.start-padFivePrime)
        end = self.end+padThreePrime
        h,s = fasta.getSequence(blastDb, self.chrom, start, end, self.strand)
        return h,s


class RichInterval:
    """An improved interval class for merging features"""
    
    attributes = ['chrom','start','end','strand','features']
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        for attrib in self.attributes:
            self.__dict__[attrib] = None
        
        if len(args)>0:
            for i,arg in enumerate(args):
                self.__dict__[self.attributes[i]] = arg
        
        if len(kw)>0:
            for attrib in kw:
                self.__dict__[attrib] = kw[attrib]
        
        if (not self.strand or self.strand=='+') and self.start>self.end:
            self.start,self.end = self.end,self.start
            self.strand = '-'
    
    @staticmethod
    def fromList(intervalList):
        """
        Static method; creates an Interval object from a list or tuple
        @param intervalList: list or tuple [start, end]
        @returns: Interval object
        """
        f = RichInterval(intervalList)
        return f
    
    @staticmethod
    def fromObject(obj):
        """
        Static method; creates an interval object from a HSP or Domain object
        
        @param obj: HSP or Domain object
        @returns: Interval object
        """
        f = RichInterval()
        f.origObject = obj
        f.start = obj.sStart
        f.end = obj.sEnd
        f.features = [obj]
        if obj.__class__==hmmer.BlockSixFrameDomain:
            f.chrom = obj.accession
            f.strand = obj.strand
        elif obj.__class__==hmmer.SixFrameDomain:
            f.chrom = obj.accession
            f.strand = obj.strand
        elif obj.__class__==blast.HSP:
            f.chrom = obj.subjectId
            f.strand = obj.strand()
        return f
    
    def toList(self):
        return self.chrom, self.start, self.end, self.strand, self.features
    
    def toDict(self):
        return dict(zip(self.attributes, self.toList()))
    
    def __repr__(self):
        return self.format % self.toDict()
    
    def overlaps(self, f):
        status = False
        if self.chrom==f.chrom and self.strand==f.strand:
            if self.start<=f.start<=self.end:
                status = True
            elif self.start<=f.end<=self.end:
                status = True
            elif f.start<=self.start<=f.end:
                status = True
        return status
    
    def strandlessOverlaps(self, f):
        status = False
        if self.chrom==f.chrom:
            if self.start<=f.start<=self.end:
                status = True
            elif self.start<=f.end<=self.end:
                status = True
            elif f.start<=self.start<=f.end:
                status = True
        return status
    
    def __cmp__(self, other):
        return cmp(self.start, other.start) or cmp(self.end, other.end)
    
    def merge(self, f):
        assert self.start<self.end
        assert f.start<f.end
        self.start = min(self.start, f.start)
        self.end = max(self.end, f.end)
        for feature in f.features:
            self.features.append(feature)


def separateBySS(features):
    featuresBySS = {}
    for f in features:
        key = (f.chrom, f.strand)
        try:
            featuresBySS[key].append(f)
        except:
            featuresBySS[key] = [f]
    
    for ss in featuresBySS:
        if ss[1]=='+':
            featuresBySS[ss] = DSU(featuresBySS[ss], lambda x: x.start)
        else:
            featuresBySS[ss] = DSU(featuresBySS[ss], lambda x: x.start)
    return featuresBySS


def merge(features):
    mergedFeatures = copy.copy(features)
    mergedFeatures.sort(key=lambda x: (x.chrom, x.strand, x.start))
    n = len(mergedFeatures)
    i = 0
    while i<n-1:
        if mergedFeatures[i].overlaps(mergedFeatures[i+1]):
            mergedFeatures[i].merge(mergedFeatures[i+1])
            del mergedFeatures[i+1]
            n -= 1
        else:
            i += 1
    return mergedFeatures


def strandlessMerge(features):
    mergedFeatures = copy.copy(features)
    mergedFeatures.sort(key=lambda x: (x.chrom, x.start))
    for f in mergedFeatures:
        f.strand = '.'
    
    n = len(mergedFeatures)
    i = 0
    while i<n-1:
        if mergedFeatures[i].strandlessOverlaps(mergedFeatures[i+1]):
            mergedFeatures[i].merge(mergedFeatures[i+1])
            del mergedFeatures[i+1]
            n -= 1
        else:
            i += 1
    return mergedFeatures


def mergeHSPs(hsps):
    features = []
    for hsp in hsps:
        features.append(RichInterval.fromObject(hsp))
    features = merge(features)
    for f in features:
        f.features.sort(key=lambda hsp: (hsp.eValue, -hsp.bitScore))
    return features


# import numpy
# def splitBitmap(interval, splits):
#     L = interval.end-interval.start+1
#     bitmap = numpy.ones(L, int)
#     for split in splits:
#         start = max(0, split.start-interval.start)
#         end = min(L-1, split.end-interval.start)
#         bitmap[start:end+1] = 0
#     intervals = []
#     for 


# class FeatureInterval:
#     """A feature interval class to be replaced by RichInterval
#     
#     The merged interval has the minimum e-value of the source intervals
#     and the maximum score. Inheritance of names also reflects the scores.
#     """
#     
#     attributes = ['name','chrom','start','end','strand','score','eValue','others']
#     format = attributesToFormat(attributes)
#     
#     def __init__(self, *args, **kw):
#         for attrib in self.attributes:
#             self.__dict__[attrib] = kw.pop(attrib, None)
#     
#     @staticmethod
#     def fromList(interval):
#         f = FeatureInterval(start=min(interval), end=max(interval))
#         return f
#     
#     @staticmethod
#     def fromObject(obj):
#         "Create a Feature object from a HSP or Domain object"
#         f = FeatureInterval()
#         f.origObject = obj
#         f.start = obj.sStart
#         f.end = obj.sEnd
#         f.eValue = obj.eValue
#         f.others = []
#         if obj.__class__==hmmer.BlockSixFrameDomain:
#             f.name = obj.domain
#             f.chrom = obj.accession
#             f.strand = obj.strand
#             f.score = obj.score
#         elif obj.__class__==blast.HSP:
#             f.name = obj.queryId
#             f.chrom = obj.subjectId
#             f.strand = obj.strand()
#             f.score = obj.bitScore
#             assert f.sStart<=f.sEnd
#         return f
#     
#     def __repr__(self):
#         return self.format % self.__dict__
#     
#     def overlaps(self, f):
#         status = False
#         if self.chrom==f.chrom and self.strand==f.strand:
#             if f.start<=self.start<=f.end:
#                 status = True
#             elif f.start<=self.end<=f.end:
#                 status = True
#             elif self.start<=f.start<=self.end:
#                 status = True
#         return status
#     
#     def merge(self, f):
#         self.start = min(self.start, f.start)
#         self.end = max(self.end, f.end)
#         self.eValue = min(self.eValue, f.eValue)
#         self.score = max(self.score, f.score)
#         
#         if self.eValue>f.eValue:
#             self.others.append(self.name)
#             self.name = f.name
#         else:
#             self.others.append(f.name)
# 
# 
# class TypedInterval(object):
#     attributes = ['chrom', 'start', 'end', 'strand', 'family']
#     format = attributesToFormat(attributes)
#     
#     def __init__(self, *args, **kw):
#         for attrib in self.attributes:
#             self.__dict__[attrib] = None
#         
#         if len(args)>0:
#             for i,arg in enumerate(args):
#                 self.__dict__[self.attributes[i]] = arg
#         
#         if len(kw)>0:
#             for attrib in kw:
#                 self.__dict__[attrib] = kw[attrib]
#         
#         if (not self.strand or self.strand=='+') and self.start>self.end:
#             self.start,self.end = self.end,self.start
#             self.strand = '-'
#     
#     @staticmethod
#     def fromList(intervalList):
#         """
#         Static method; creates an Interval object from a list or tuple
#         @param intervalList: list or tuple [start, end]
#         @returns: Interval object
#         """
#         f = Interval(intervalList)
#         return f
#     
#     def toList(self):
#         return self.chrom, self.start, self.end, self.strand, self.family
#     
#     def toDict(self):
#         return dict(zip(self.attributes, self.toList()))
#     
#     def __repr__(self):
#         return self.format % self.toDict()
#     
#     def overlaps(self, f):
#         status = False
#         if self.strand==f.strand and self.family==f.family:
#             if self.start<=f.start<=self.end:
#                 status = True
#             elif self.start<=f.end<=self.end:
#                 status = True
#             elif f.start<=self.start<=f.end:
#                 status = True
#         return status
#     
#     def merge(self, f):
#         self.start = min(self.start, f.start)
#         self.end = max(self.end, f.end)
