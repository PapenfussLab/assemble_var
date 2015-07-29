"""
genscan module
"""

import fasta
from useful import smartopen
from mungoCore import *




class Predicted(AbstractFeature):
    """Predicted gene class"""
    
    attributes = [
        'gene_exon','type','strand','start','end','length',
        'Fr','Ph','I_Ac','Do_T','CodRg','P','score']
    converters = []
    format = attributesToFormat(attributes)
    
    def __init__(self, *args, **kw):
        super(Predicted, self).__init__(*args,**kw)
    
    def convert(self):
        intVars = ['start','end','length','Fr','Ph','I_Ac','Do_T','CodRg']
        floatVars = ['score', 'P']
        casts = [int]*len(intVars) + [float]*len(floatVars)
        for item,cast in zip(intVars+floatVars, casts):
            try:
                self.__dict__[item] = cast(self.__dict__[item])
            except:
                self.__dict__[item] = None
    
    def __repr__(self):
        d = {}
        for k,v in self.__dict__.iteritems():
            d[k] = v
        return self.format % d


def load(iFileHandle):
    """Load genscan predictions.
    
    Arguments:
    iFileHandle -- Input file or filename.
    
    Return values:
    data -- Annotation data (a list of lists, each list in one gene)
    
    """
    return load_full(iFileHandle)[0]


def load_preprocessed(iFileHandle):
    """Load genscan predictions when predictions have been preprocessed
    and only contain the gene prediction lines
    
    Arguments:
    iFileHandle -- Input file or filename.
    
    Return values:
    data -- Annotation data (a list of lists, each list in one gene)
    
    """
    iFile = smartopen(iFileHandle)
    data = {}
    
    skipState = 'Slice no. '
    
    state = None
    for line in iFile:
        line = line.strip()
        if line:
            if skipState in line:
                pass
            else:
                tokens = line.split()
                d = Predicted(tokens)
                gene = int(d.gene_exon.split('.')[0])
                try:
                    data[gene].append(d)
                except KeyError:
                    data[gene] = [d]
    
    data = data.items()
    data.sort()
    data = [x[1] for x in data]
    
    return data
    

def load_full(iFileHandle):
    """Load genscan predictions.
    
    Arguments:
    iFileHandle -- Input file or filename.
    
    Return values:
    data -- Annotation data (a list of lists, each list in one gene)
    proteins -- Predicted proteins (a list of tuples (header, sequence))
    meta -- Meta-data in first 8 lines of genscan output
    
    """
    iFile = smartopen(iFileHandle)
    data = {}
    proteins = []
    meta = []
    
    startPredState = '----- ---- - ------ ------ ---- -- -- ---- ---- ----- ----- ------'
    endPredState = 'Predicted peptide sequence(s):'
    skipState = 'Slice no. '
    metaState = 'GENSCAN 1.0'
    
    state = None
    for line in iFile:
        line = line.strip()
        
        if metaState in line:
            state = 'meta'
        if line==startPredState:
            state = 'pred'
        elif line=='NO EXONS/GENES PREDICTED IN SEQUENCE':
            state = 'fail'
        elif line==endPredState:
            state = 'prot'
        elif skipState in line:
            state = 'skip'
        else:
            if state=='meta':
                if line:
                    meta.append(line)
            elif state=='pred':
                if line:
                    tokens = line.split()
                    d = Predicted(tokens)
                    gene = int(d.gene_exon.split('.')[0])
                    try:
                        data[gene].append(d)
                    except KeyError:
                        data[gene] = [d]
            elif state=='prot':
                break
            elif state=='fail':
                return [], [], ''
    
    if state=='prot':
        proteins = fasta.load_mfa(iFile)
    
    data = data.items()
    data.sort()
    data = [x[1] for x in data]
    
    return data, proteins, meta



def tokens2dict(tokens):
    """Convert a genscan token list to a dictionary."""
    if tokens[1] in ['Prom','PlyA']:
        d = dict(zip(endFields, tokens))
        for attrib in Predicted.attributes:
            if not d.has_key(attrib):
                d[attrib] = '.'
    else:
        d = dict(zip(Predicted.attributes, tokens))
        
    if d['strand']=='-':
        d['start'],d['end'] = d['end'],d['start']
    
    converters = {'start': int, 'end': int, 'score': float}
    for attrib,convert in converters.iteritems():
        try:
            d[attrib] = convert(d[attrib])
        except KeyError:
            print tokens
    
    return d
