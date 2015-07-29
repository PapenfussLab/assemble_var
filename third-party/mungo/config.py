"""
config class & objects
"""

import os


class Config:
    """Config class"""
    def __init__(self, *args, **kw):
        """Constructor
        
        Arguments:
        name: Name (default: '')
        folder: Root directory for species (say) (default: '')
        """
        if len(args)>0:
            self.name = args[0]
        else:
            self.name = ''
        self.folder = kw.pop('folder', '')
        for attr in kw:
            self.__dict__[attr] = kw[attr]
    
    def add(self, attr, value, relative=True):
        if relative:
            self.__dict__[attr] = os.path.join(self.folder, value)
        else:
            self.__dict__[attr] = value
    
    def __repr__(self):
        return str(self.__dict__)


homeDir = os.environ['HOME']
dataDir = os.path.join(homeDir, 'databases')


# RefSeq
refseqData = Config('RefSeqData', folder=os.path.join(dataDir, 'refseq'))
refseqData.add('accession2geneid', 'accession2geneid')
refseqData.add('all', 'blastdb/all')

# refseqData.mammal = Config('Mammalian RefSeq', folder=refseqData.folder)
# refseqData.mammal.add('rnaBlastDb', 'blastdb/vertebrate_mammalian.rna')
# refseqData.mammal.add('proteinBlastDb', 'blastdb/vertebrate_mammalian.pep')
# 
# refseqData.protozoa = Config('Protozoa RefSeq', folder=refseqData.folder)
# refseqData.protozoa.add('proteinBlastDb', 'blastdb/protozoa.pep')


# Opossum
opossum = Config('Opossum', folder=os.path.join(dataDir, 'opossum'))
opossum.add('genomeBlastDb', 'assembly/blastdb/assembly')
opossum.add('chrSizes', 'assembly/chrSizes.txt')
opossum.add('chrSizesWoChr', 'assembly/chrSizesWoChr.txt')


# Platypus
platypus = Config('Platypus', folder=os.path.join(dataDir, 'platypus'))
platypus.add('genomeBlastDb', 'assembly/assembly')
platypus.add('chrSizes', 'assembly/chrSizes.txt')

platypus.ensembl = Config('Platypus ensembl predictions', 
    folder=os.path.join(platypus.folder, 'ensembl/latest'))
platypus.ensembl.add('peptides', 'Ornithorhynchus_anatinus.OANA5.46.pep.all.fa')
# platypus.ensembl.add('cdna', 'Ornithorhynchus_anatinus.OANA5.46.cdna.all.fa')


# Human
human = Config('Human', folder=os.path.join(dataDir, 'human'))
human.add('genomeBlastDb', 'blastdb/hg17')
human.add('chrSizes', 'chrom_sizes.txt')

# Human RefSeq
human.refseq = Config('Human RefSeq',
    folder=os.path.join(human.folder, 'refseq/latest'))
human.refseq.add('flat', 'refFlat.txt')
human.refseq.add('link', 'refLink.txt')
human.refseq.add('gene', 'refGene.txt')
human.refseq.add('proteinBlastDb', 'refseq')


# Mouse
mouse = Config('Mouse', folder=os.path.join(dataDir, 'mouse'))

# Mouse RefSeq
mouse.refseq = Config('Mouse RefSeq', 
    folder=os.path.join(mouse.folder, 'refseq/v19'))
mouse.refseq.add('flat', 'refFlat.txt')
mouse.refseq.add('link', 'refLink.txt')
mouse.refseq.add('gene', 'refGene.txt')
mouse.refseq.add('proteinBlastDb', 'refseq')


# Parasites
Pf = Config('Plasmodium falciparum',
    folder=os.path.join(dataDir, 'apicomplexa/plasmodium/falciparum'))
Pf.add('proteinBlastDb', 'proteins')

Py = Config('Plasmodium yoelli',
    folder=os.path.join(dataDir, 'apicomplexa/plasmodium/yoelli'))
Py.add('proteinBlastDb', 'proteins')

Lmajor = Config('Leishmania major', 
    folder=os.path.join(dataDir, 'leishmania'))
Lmajor.add('proteinBlastDb', 'proteins')
