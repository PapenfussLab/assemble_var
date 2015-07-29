"""
Sequence module

Useful sequence analysis functions
"""

import random, re, copy, math
import hmmer


nucleotides = ['A','T','G','C']

aminoAcids = ['A','C','E','D','G','F','I','H','K','M','L','N','Q','P','S','R','T','W','V','Y']

baseComplement = {
    'A':'T', 'C':'G', 'G':'C', 'T':'A',
    'a':'t', 'c':'g', 'g':'c', 't':'a',
    'N':'N', 'n':'n', 'X':'X', '-':'-',
    'R':'Y', 'Y':'R', 'M':'K', 'K':'M',
    'r':'y', 'y':'r', 'm':'k', 'k':'m',
    'W':'W', 'S':'S', 'B':'V', 'V':'B',
    'w':'w', 's':'s', 'b':'v', 'v':'b',
    'D':'H', 'H':'D', 'd':'h', 'h':'d'}

standardCodonTable = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'NNN': 'X',                 
    'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L',
    'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
    'tat': 'Y', 'tac': 'Y', 'taa': '*', 'tag': '*',
    'tgt': 'C', 'tgc': 'C', 'tga': '*', 'tgg': 'W',
    'ctt': 'L', 'ctc': 'L', 'cta': 'L', 'ctg': 'L',
    'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P',
    'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q',
    'cgt': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R',
    'att': 'I', 'atc': 'I', 'ata': 'I', 'atg': 'M',
    'act': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T',
    'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K',
    'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R',
    'gtt': 'V', 'gtc': 'V', 'gta': 'V', 'gtg': 'V',
    'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A',
    'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E',
    'ggt': 'G', 'ggc': 'G', 'gga': 'G', 'ggg': 'G'
}

uniformNucleotideDist = {'G': 0.25, 'C': 0.25, 'A': 0.25, 'T': 0.25}

toyNucleotideDist = {'G': 0.37, 'C': 0.37, 'A': 0.5*(1-2*0.37), 'T': 0.5*(1-2*0.37)}

aminoAcidDist = {
    'A': 0.075520,
    'C': 0.016973,
    'D': 0.053029,
    'E': 0.063204,
    'F': 0.040762,
    'G': 0.068448,
    'H': 0.022406,
    'I': 0.057284,
    'K': 0.059398,
    'L': 0.093399,
    'M': 0.023569,
    'N': 0.045293,
    'P': 0.049262,
    'Q': 0.040231,
    'R': 0.051573,
    'S': 0.072214,
    'T': 0.057454,
    'V': 0.065252,
    'W': 0.012513,
    'Y': 0.031985,
}

lnAminoAcidDist = {}
for aa in aminoAcidDist:
    lnAminoAcidDist[aa] = math.log10(aminoAcidDist[aa])


def uniformRandomLetter(alphabet=nucleotides):
    return random.choice(alphabet)


def uniformRandomSequence(L, alphabet=nucleotides):
    seq = []
    for i in xrange(L):
        seq.append(uniformRandomLetter(alphabet))
    return ''.join(seq)


def randomLetter(distribution=toyNucleotideDist):
    r = random.random()
    cumsum = 0.0
    for letter,freq in toyNucleotideDist.iteritems():
        cumsum += freq
        if r<=cumsum:
            return letter


def randomSequence(L, distribution=toyNucleotideDist):
    seq = []
    for i in xrange(L):
        seq.append(randomLetter(distribution))
    return ''.join(seq)


def transcribe(dna):
    """Return DNA string as RNA string."""
    return dna.replace('T', 'U').replace('t','u')


def reverse(s):
    """Return the sequence string in reverse order."""
    letters = list(s)
    letters.reverse()
    return ''.join(letters)


def complement(s):
    """Return the complementary sequence string."""
    letters = list(s)
    letters = [baseComplement[base] for base in letters]
    return ''.join(letters)


def reverseComplement(s):
    """Return the reverse complement of the sequence string."""
    letters = list(s)
    letters = [baseComplement[base] for base in letters]
    letters.reverse()
    return ''.join(letters)



def AT(s):
    "Return the percentage of dna composed of A+T"
    s = s.upper()
    n = s.count("A") + s.count("T")
    nn = s.count("N")
    denom = len(s)-nn
    if denom==0:
        return 0.0
    else:
        return 100.0*n / denom


def GC(s):
    "Return the percentage of dna composed of G+C."
    s = s.upper()
    n = s.count("G") + s.count("C")
    nn = s.count("N")
    denom = len(s)-nn
    if denom==0:
        return 0.0
    else:
        return 100.0*n / denom


class FrameError(Exception):
    message = 'abs(frame)>3!!!'


def codons(s, remainder=False):
    "Return list of codons for the dna string & optionally return the remainder."
    end = len(s) - len(s) % 3
    codons = [s[i:i+3] for i in range(0, end, 3)]
    if remainder:
        rem = s[end:]
        return codons, rem
    else:
        return codons


def codonIterator(s, remainder=False):
    "Returns an iterator that yields codons & optionally the remainder."
    L = len(s)
    i = 0
    while True:
        if i+3<=L:
            yield s[i:i+3]
        elif L-i>0:
            break
        elif L-i==0:
            return
        i += 3
    
    if remainder:
        yield s[i:]


def translate(s, codonTable=standardCodonTable, remainder=False, frame=0):
    s2,rem = codons(s[frame:].upper(), remainder=True)
    protein = []
    for codon in s2:
        try:
            aa = codonTable[codon]
        except KeyError:
            aa = 'X'
        protein.append(aa)
    if remainder:
        return ''.join(protein), rem
    else:
        return ''.join(protein)


def sixFrameTranslationIter(s, codonTable=standardCodonTable):
    """Returns a 6 frame translation iterator.
    
    @param s: DNA sequence
    @param codonTable: Codon translation table (default=standardCodonTable)
    @return: Iterator yielding (frame, translation)
    """
    # First 3 frames
    for frame in xrange(3):
        yield frame+1, translate(s[frame:])
    
    # Reverse strand
    s2 = reverseComplement(s)
    for frame in xrange(3):
        yield -frame-1, translate(s2[frame:])


def sixFrameTranslation(seq, codonTable=standardCodonTable):
    """Returns a dictionary containing a 6 frame translation.
    @param seq: DNA sequence
    @param codonTable: Codon translation table (default=standardCodonTable)
    """
    trans = {}
    for frame,pep in sixFrameTranslationIter(seq):
        trans[frame] = pep
    return trans


def extractOrfsIter(seq, minLen=20, pattern='\*|X{200,}', verbose=False):
    """Returns an ORF extracting iterator
    @param seq: DNA sequence
    @param minLen: Minimum ORF length (default=20)
    @param pattern: Stop and break pattern (default="*|X{200,}")
    @return: Iterator yielding (i,gStart,gEnd,ORF)
    """
    L = len(seq)
    regex = re.compile(pattern)
    
    i = 0
    sixFrameIter = sixFrameTranslationIter(seq)
    for frame,p in sixFrameIter:
        if verbose: print frame
        matchIter = regex.finditer(p)
        
        # As though there is a stop to the left of 0
        start = -1
        for match in matchIter:
            end = match.start()
            orf = p[start+1:end]
            if len(orf)>=minLen:
                i += 1
                gStart,gEnd = hmmer.convertSixFrameToGenomic(start+2, end, frame, L)
                yield i, gStart, gEnd, orf
            start = copy.copy(end)
        
        # As though there is a stop to the right of len(p)-1
        end = len(p)
        orf = p[start+1:end]
        if len(orf)>=minLen:
            i += 1
            gStart,gEnd = hmmer.convertSixFrameToGenomic(start+2, end, frame, L)
            yield i, gStart, gEnd, orf


def extractOrfs(seq, baseHeader='', minLen=20, pattern='\*|X{200,}'):
    """Returns all ORFs >minLen from a sequence.
    
    @param seq: DNA sequence
    @param baseHeader: Base name for header (default='')
    @param minLen: Minimum ORF length (default=20)
    @param pattern: Stop and break pattern (default="*|X{200,}")
    @return: List of [header, ORF]
    """
    if baseHeader: baseHeader = baseHeader + '.'
    data = []
    orfIter = extractOrfsIter(seq, minLen=minLen, pattern=pattern)
    for i,gStart,gEnd,orf in orfIter:
        h = '%s%i %i-%i Length=%i' % (baseHeader,i,gStart,gEnd,len(orf))
        data.append([h, orf])
    return data


def calcFrequency(seq, alphabet=aminoAcids):
    freq = {}
    for symbol in alphabet: freq[symbol] = 0
    for symbol in seq: 
        _symbol = symbol.upper()
        try:
            freq[_symbol] += 1
        except:
            pass
    
    total = float(sum(freq.values()))
    for symbol in alphabet: freq[symbol] /= total
    return freq

# ----------------------------------------------------------------------

from useful import warnDeprecated


def reverse_complement(s):
    """@deprecated: use sequence.reverseComplement"""
    warnDeprecated('sequence.reverse_complement, use reverseComplement.')
    return reverseComplement(s)


