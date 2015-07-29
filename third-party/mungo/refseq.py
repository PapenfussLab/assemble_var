"""
RefSeq module
"""

import os, sys
import httplib, ftplib
from config import human


def updateFromUCSC(build='hg18', oDir='.'):
    filenames = [
        #'refGene.txt.gz',
        #'refLink.txt.gz',
        #'refFlat.txt.gz',
        'refGene.sql',
        'refLink.sql',
        'refFlat.sql',
    ]
    
    conn = httplib.HTTPConnection('hgdownload.cse.ucsc.edu')
    for filename in filenames:
        print filename
        url = '/goldenPath/%s/database/%s' % (build, filename)
        conn.request('GET', url)
        resp = conn.getresponse()
        if resp.status==200:
            oFilename = '%s/%s' % (oDir, filename)
            f = open(oFilename, 'w')
            f.write(resp.read())
            f.close()
        else:
            print 'Error downloading', url, resp.status, resp.reason
    conn.close()


def updateFromNCBI(spp='human'):
    data = {
        'human': {
            'iFilename': '/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz',
            'oFilename': './human.proteins.faa.gz'
        }
    }
    
    ftp = ftplib.FTP('ftp.ncbi.nih.gov')
    ftp.login()            
    f = open(data[spp]['oFilename'], 'w')
    ftp.retrbinary('RETR %s' % data[spp]['iFilename'], f.write)
    ftp.quit()
    f.close()


def constructDb(iDir='.'):
    pass


if __name__=='__main__':
    updateFromUCSC()
    updateFromNCBI()
    