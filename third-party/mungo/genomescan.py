"""
genomescan module
"""

import os, sys
import time
import fasta, genscan, sequence
from useful import smartopen
from useful import multipart


def genomescan(dna, protein, description='', email='', oFileHandle=None, 
  proxy='wehiproxy.alpha.wehi.edu.au', proxyPort=3128, debug=False):
    """
    
    @param dna: DNA sequence in fasta format.
    @param protein: protein sequence in fasta format.
    @param description: Description of job (default: '').
    @param email: Email address (default: '').
    @param oFileHandle: Output file or filename (default: None/no output)
    @param proxy: Proxy setting (default: WEHI proxy)
    @param proxyPort: Proxy port setting (default: WEHI proxy port)
    @param debug: Debug output (default: False)
    @returns: HTML output of GenomeScan webserver.
    """
    genomeScanURL = 'http://genes.mit.edu/cgi-bin/genomescanw.cgi'
    
    if debug:
        print dna
        print protein
    
    # Do genomescan prediction using webserver
    fields = [
        ('-g', protein), 
        ('-s', dna),
        ('-o', 'Vertebrate'),
        ('-n', description),
        ('-p', 'Predicted peptides only'),
        ('-a', email),
    ]
    files = [('-u', '', ''), ('-v', '', '')]
    
    try:
        html = multipart.post(genomeScanURL, fields, files, proxy, proxyPort)
    except multipart.FormSubmissionException:
        print >> sys.stderr, "*** %s submission failed. Retry later" % description
        return
    except Exception, e:
        print e
        sys.exit('Argh!')
    
    if oFileHandle:
        oFile = smartopen(oFileHandle, 'w')
        print >> oFile, html
        oFile.close()
    
    return html


def extractSeq(feature, blastDb, dx, dy):
    """Extract the translated sequence of a feature and DNA 
    sequence of the surrounding the genomic region.
    
    @param feature: Feature object. Mandatory attributes: accession, sStart, sEnd.
    @param blastDb: Blast database.
    @param dx: Length of sequence to extract upstream.
    @param dy: Length of sequence to extract downstream.
    @returns: tuple of fasta strings (DNA, protein).
    """
    # Extract hmmer hit & translate
    header,seq = fasta.getSequence(blastDb, feature.accession, 
        start=feature.sStart, end=feature.sEnd, strand=feature.strand)
    protein = '\n'.join(['>' + header, sequence.translate(seq)])
    
    if feature.strand=='+':
        start = feature.sStart-dx
        end = feature.sEnd+dy
    else:
        start = feature.sStart-dy
        end = feature.sEnd+dx
    
    if start<0 or end>5000000:
        raise Exception('Out of block bounds.')
    
    # Extract surrounding DNA sequence
    header,seq = fasta.getSequence(blastDb, feature.accession, 
        start=start, end=end)
    header = '%(accession)s:%(sStart)s-%(sEnd)s' % feature.__dict__
    dna = '\n'.join(['>' + header, seq])
    
    return dna, protein


def genomescanFromFeature(feature, blastDb, dx=15000, dy=7000, 
  oFileHandle=None, debug=False):
    """Predict a gene in the neighbourhood of a known feature.
    
    @param feature: Feature object. Mandatory attributes: accession, sStart, sEnd.
    @param blastDb: Blast database.
    @param dx: Length of sequence to extract upstream (default: 20000)
    @param dy: Length of sequence to extract downstream (default: 7000)
    @param oFileHandle: Output file or filename (default: None/no output)
    @param debug: Debug output (default: False)
    @returns: HTML output of GenomeScan webserver.
    """
    description = feature.domain
    if debug: print '%s  L=%i  time="%s"' \
      % (description, feature.sEnd-feature.sStart+1+dx+dy, time.ctime())
    dna,protein = extractSeq(feature, blastDb, dx, dy)
    html = genomescan(dna, protein, description=description, \
        oFileHandle=oFileHandle, debug=debug)
    return html


def parseGenomeScanOutput(html, debug=False):
    """A simple stateful parser for genomescan output.
    
    @param html: Raw HTML output of genomescan webserver.
    @param debug: Debug output (default: False)
    @returns: Tuple (annotation, predicted proteins in fasta format).
    """
    class Start:
        pass
    
    class AnnotationStart:
        token = '<b>Predicted genes/exons:'
    
    class AnnotationEnd:
        token = '<a href=#exp2>explanation</a>'
    
    class NoAnnotation:
        token = '</b><b>NO EXONS/GENES PREDICTED IN SEQUENCE'
    
    class PeptideStart:
        token = '<b>Predicted peptide sequence(s):'
    
    class PeptideEnd:
        token = '<a name=exp1>'
    
    annot = []
    header = ''
    seq = []
    
    state = Start
    htmlIter = iter(html)
    for line in htmlIter:
        line = line.strip()
        
        if line==AnnotationStart.token:
            state = AnnotationStart
            continue
        elif line==AnnotationEnd.token:
            state = AnnotationEnd
            continue
        elif line==NoAnnotation.token:
            state==NoAnnotation
            break
        elif line==PeptideStart.token:
            state = PeptideStart
            htmlIter.next()
            htmlIter.next()
            continue
        elif line==PeptideEnd.token:
            state = PeptideEnd
            break
        
        if state==AnnotationStart:
            if line[0:7]=='</b><b>':
                annot.append(line[7:])
            if debug: print annot
        elif state==PeptideStart:
            seq.append(line)
            if debug: print seq
    
    return '\n'.join(annot), '\n'.join(seq)


def parseGenscan(lines, skipLines=2):
    """Parse lines in genscan or genomescan output.
    
    @param lines: Lines from genscan/genomescan annotation.
    @param skipLines: Number of lines to skip at start (default 3).
    """
    data = []
    for line in lines[skipLines:]:
        line = line.strip()
        if line:
            tokens = line.split()
            d = genscan.Predicted(tokens)
            # gene = int(d.gene_exon.split('.')[0])
            try:
                data[-1].append(d)
            except:
                data.append([d])
    return data
