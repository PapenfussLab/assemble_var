# Pipeline for the assembly of VAR genes

This pipeline aims to assemble Plasmodium VAR genes from RNA sequence. It uses a number of third party programs including Trim Galore, Subread, Pear, Khmer, Oases and Cap3.  Here I give a brief outline of how to use it.

The main program is **assembly_var.py** it takes as input:
* A reference fasta file with all the known contaminant genomes. For instance I used *human*, *Plasmodium vivax* and *Plasmodium falciparum*. 
* A number of gff files indicating which regions of the reference are of interest. I used the VAR gene regions of *Plasmodium vivax* and *Plasmodium falciparum* as found in PlasmoDB.
* Two paired end read files in fastq format

It ouputs:
* A fasta file with contigs. This is stored in the output directory in a file named **transcripts_assembled_final.fa**.  Further filtering may then be required to remove any remaining contaminants.


### assemble_var.py
```

Options:
  -h, --help            show this help message and exit
  -r READ1, --read1=READ1
                        first set of read pairs
  -R READ2, --read2=READ2
                        second set of read pairs
  --reference=REFERENCE
                        a fasta file containing the reference genomes we want
                        to filter out. If using multiple different references
                        they will need to be combined into one fasta file. i.e
                        using cat
  --index=REF_INDEX     the location of the index files. Usefule if performing
                        multiple assemblies as this does not need to be
                        recomputed each time.
  -v VAR_FILES, --var=VAR_FILES
                        var files that contain var regions we want to keep -
                        that is a subset of reference
  -o OUTPUTDIR, --outputdir=OUTPUTDIR
                        output directory. This will be created if it doesn't
                        already exist.
  --ins_length=INS_LENGTH
                        the insert length passed to oases
  --adapter=ADAPTER     the adapter sequence to be passed to Trim Galore.
                        Defaults to trim galore's default if not supplied.
  --verbose             turns on more detailed output
  --pear                merge read pairs that overlap before assembly.
  --norm                perform digital normalisation which decreases the
                        computational time required for assembly.
  --soap                assemble reads using soapDenovo-trans.

```
#### Usage Example

```python
python assemble_var.py -r /path/to/read1.fastq -R /path/to/read2.fastq --reference /path/to/reference.fasta -v /path/to/Pfalciparum_VAR.gff -o /path/to/output_directory/ --norm
```
```
