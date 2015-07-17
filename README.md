# Pipeline for the assembly of VAR genes

This pipeline aims to assemble Plasmodium VAR genes from RNA sequence. It uses a number of third party programs including Trim Galore, Subread, Pear, Khmer, Oases and Cap3.  Here I give a brief outline of how to use it.

The main program is **assembly_var.py** it takes as input:
* A reference fasta file with all the known contaminant genomes. For instance I used *human*, *Plasmodium vivax* and *Plasmodium falciparum*. 
* A number of gff files indicating which regions of the reference are of interest. I used the VAR gene regions of *Plasmodium vivax* and *Plasmodium falciparum* as found in PlasmoDB.
* Two paired end read files in fastq format

It ouputs:
* A fasta file with contigs. Further filtering may then be required to remove any remaining contaminants.


### assemble_var.py
```python

    parser.add_option("-r", "--read1", dest="read1",
        help="first set of read pairs")

    parser.add_option("-R", "--read2", dest="read2",
        help="second set of read pairs")

    parser.add_option("", "--reference", dest="reference"
        , default=None
        , help=("a fasta file containing the reference"
            + " genomes we want to filter out"))

    parser.add_option("", "--index", dest="ref_index", default=False,
        help=("the location of the index files."
            + " Usefule if performing multiple assemblies as this"
            + " does not need to be recomputed each time."))

    parser.add_option('-v', '--var',
                  type='string', action='append',
                  default=[], dest='var_files',
                  help=('var files that contain var regions we want to keep'
                    + ' - that is a subset of reference'))

    parser.add_option("-o", "--outputdir", dest="outputdir",
        help="output directory")

    parser.add_option("", "--ins_length", dest="ins_length"
        , default=False, help="the insert length passed to oases")

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("", "--pear", action="store_true", dest="pear"
        , default=False, help="merge read pairs that overlap before oases.")

    parser.add_option("", "--norm", action="store_true", dest="norm"
        , default=False, help="merge read pairs that overlap before oases.")
```
###### -r --read1
The first fastq file
###### -R --read2
The second fastq file
###### ssg
