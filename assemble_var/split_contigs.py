import sys, os
from mungo.fasta import FastaReader
from collections import defaultdict
from optparse import OptionParser
from subprocess import check_call

from third_party_runners import run_blast

BLAST_NT_DB = "/home/users/allstaff/tonkin-hill.g/find_var_genes/assemble_var/data/nt/nt"

def reNameContigs(contig_file, fileName, outputdir):
    renamed = outputdir + fileName + "renamed.fa"
    with open(renamed , 'w') as outfile:
        for h,s in FastaReader(contig_file):
            h = h.split()[0]
            if len(h.split("_")) > 4: #oases transcript
                h = h.split("_")
                h = "_".join([h[0], h[1], h[3]])
            outfile.write(">" + h + "\n")
            outfile.write(s + "\n")
    return renamed

def filter_length(contig_file, length_filter, fileName, outdir, verbose):
  length_file = outdir + fileName + "lenFilt.fa"

  short_count = 0
  with open(length_file,'w') as outfile:
    for h,s in FastaReader(contig_file):
      if len(s)<length_filter:
        short_count+=1
      else:
        outfile.write(">" + h + "\n")
        outfile.write(s + "\n")

  if verbose:
    print short_count, " contigs removed as too short..."

  return length_file


def get_contaminants(fasta_ref_files, contig_file, fileName, percent_overlap
    , outdir, verbose):
  #first get list of contigs
  contigs = {}
  for h,s in FastaReader(contig_file):
      contigs[h] = s

  if verbose:
    print ("Number of contigs before contaminant filtering: "
      , len(contigs.keys()))

  #now run blast against the reference files which we want not to be
  #present in the data i.e. human
  blast_files = []
  for reference in fasta_ref_files:
      blast_files.append(run_blast(reference, contig_file, outdir
          , verbose))

  #now iterate through blast results file removing contigs that have to
  #high a proportion of hits
  bad_contigs = set()
  for blast_file in blast_files:
      blast_name = os.path.splitext(os.path.basename(blast_file))[0]
      with open(blast_file, 'r') as bfile:
          for line in bfile:
              tokens = line.strip().split()
              name = tokens[0]
              overlap = int(tokens[3])/float(len(contigs[name]))
              if overlap > percent_overlap:
                  #we don't want this contig
                  print "removing", name, "overlapped", blast_name
                  bad_contigs.add(name)

  #now write out a fasta file of contaminant sequences
  contaminant_file = outdir + fileName + "contaminants.fa"
  with open(contaminant_file, 'w') as outfile:
    for contig in bad_contigs:
        outfile.write(">" + contig + "\n")
        outfile.write(contigs[contig] + "\n")

  #now write contigs without contaminants to a file
  non_contaminant_file = outdir + fileName + "Non_contaminants.fa"
  with open(non_contaminant_file, 'w') as outfile:
      for contig in contigs:
        if contig not in bad_contigs:
          outfile.write(">" + contig + "\n")
          outfile.write(contigs[contig] + "\n")

  if verbose:
    print ("Number of contigs after filtering: "
      , len(contigs.keys())-len(bad_contigs))

  return non_contaminant_file, contaminant_file


def get_rask_var(raskFasta, contig_file, fileName, outdir, verbose):

  #run blast against the rask VAR genes
  blast_out = run_blast(raskFasta, contig_file, outdir, verbose)

  rask_hits = set()
  with open(blast_out , 'r') as blastfile:
      for line in blastfile:
          rask_hits.add(line.split()[0])

  rask = outdir + fileName + "_rask.fa"
  non_rask = outdir + fileName + "_nonrask.fa"

  count_rask = 0
  count_non_rask = 0

  with open(rask, 'w') as raskout:
      with open (non_rask, 'w') as nonraskout:
          for h,s in FastaReader(contig_file):
              if h in rask_hits:
                  raskout.write(">"+h+"\n")
                  raskout.write(s+"\n")
                  count_rask += 1
              else:
                  nonraskout.write(">"+h+"\n")
                  nonraskout.write(s+"\n")
                  count_non_rask +=1

  if verbose:
    print count_rask, " annotated to rask DB..."
    print count_non_rask, " remaining"

  return rask, non_rask

def annotate_w_ntDB(contig_file, fileName, outdir, verbose):

  blastOut = outdir + fileName + "nonRaskBlast.txt"


  #first run a special blast using the nt database
  blast_cmd = ("blastn "
      + "-evalue 10 "
      + """-outfmt "6 qseqid sseqid  stitle length pident qstart qend sstart send evalue" """
      # + "-num_alignments " + str(num_hits) + " "
      + "-num_threads 10 -max_target_seqs 3 "
      + "-db " + BLAST_NT_DB + " "
      + "-query " + contig_file + " "
      + "-out " + blastOut
      )
  if verbose:
      print blast_cmd

  check_call(blast_cmd, shell=True)

  #now retrieve annotation information
  contigs = defaultdict(str)
  contigs_perID = defaultdict(float)
  with open(blastOut, 'r') as blastfile:
    for line in blastfile:
      line = line.strip().split("\t")
      contigs[line[0]] = (contigs[line[0]]
        +  " [" + line[2]
        + "_alignLen_" + line[3]
        + "_perID_" + line[4]
        + "] ")
      contigs_perID[line[0]] = max(contigs_perID[line[0]], float(line[4]))

  #now re-write the fasta file with the annotations in the headers
  annotated = outdir + fileName + "nonRask_annotated.fa"
  unknown_blastOut = outdir + fileName + "ForManualInspection.fa"

  with open(annotated, 'w') as outfileKnown:
    with open(unknown_blastOut, 'w') as outfileUnknoen:
      for h,s in FastaReader(contig_file):
        if h in contigs:
          if contigs_perID[h] > 0.97:
            outfileKnown.write(">" + h + " " + contigs[h] + "\n")
            outfileKnown.write(s + "\n")
            continue
        if h in contigs:
          outfileUnknoen.write(">" + h + " " + contigs[h] + "\n")
        else:
          outfileUnknoen.write(">" + h + " none\n")
        outfileUnknoen.write(s + "\n")


  return annotated


def split_contigs(contig_file, contaminant_refs, raskDB, min_length
  , perThreshold, outdir, verbose):

  #get contig filenames for easy naming of new files
  fileName = os.path.splitext(os.path.basename(contig_file))[0]

  renamed = reNameContigs(contig_file, fileName, outdir)

  filtered = filter_length(renamed, min_length, fileName, outdir, verbose)

  non_contaminant_file, contaminant_file = get_contaminants(contaminant_refs
    , filtered, fileName, perThreshold
    , outdir, verbose)

  rask, non_rask = get_rask_var(raskDB, non_contaminant_file, fileName, outdir, verbose)

  annotated = annotate_w_ntDB(non_rask, fileName, outdir, verbose)

  return


def main():
    parser = OptionParser()

    parser.add_option("", '--ref',
                  type='string', action='append',
                  default=[], dest='ref_files',
                  help=('reference files that contain the sequence we want to '
                    + 'filter the contigs with'))

    parser.add_option("-o", "--outputdir", dest="outputdir",
        help="output directory")

    parser.add_option("-p", "--perThreshold", dest="percent", type=float
        , help=("percent of contig length that is allowed to overlap with "
            +"80%->0.8 reference"))

    parser.add_option("-c", "--contig", dest="contig_file",
        help="fas file of contigs that are to be filtered")

    parser.add_option("", "--min_length", dest="min_length", type=int
        , help="the minimum length a contig must be in order to keep it")

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("", "--raskDB", dest="raskDB",
        help="fasta file containing the VAR gene sequences of Rask et al")


    (options, args) = parser.parse_args()

    split_contigs(options.contig_file, options.ref_files, options.raskDB
        , options.min_length, options.percent
        , options.outputdir, options.verbose)


if __name__ == '__main__':
    main()



















