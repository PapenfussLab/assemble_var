# Pipeline for the assembly of VAR genes

This pipeline aims to assemble Plasmodium VAR genes from RNA sequence. It uses a number of third party programs including Trim Galore, Subread, Pear, Khmer, Oases and Cap3.  Here I give a brief outline of how to use it.

The main program is **assembly_var.py** it takes as input:
* A reference fasta file with all the known contaminant genomes. For instance I used *human*, *Plasmodium vivax* and *Plasmodium falciparum*. 
* A number of gff files indicating which regions of the reference are of interest. I used the VAR gene regions of *Plasmodium vivax* and *Plasmodium falciparum* as found in PlasmoDB.
* Two paired end read files in fastq format

It ouputs:
* A fasta file with contigs. Further filtering may then be required to remove any remaining contaminants.

## assemble_var.py

### test
