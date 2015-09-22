from mungo.sequence import sixFrameTranslation
from mungo.fasta import FastaReader
import os, sys
from optparse import OptionParser
from subprocess import check_call

HMMER3 = "/home/users/allstaff/tonkin-hill.g/rask_based_block_finder2.0/third-party/hmmer-3.1b1/src/hmmsearch"

def split_easy_from_hard(inputfile, outputdir):
  seqCount = 0
  badSeqs = 0
  bad_lengths = []

  output_file = (outputdir + os.path.splitext(os.path.basename(inputfile))[0]
    + "_translatedL2stops.fa")

  with open(output_file + "_BadSeqs", 'w') as badfile:
    with open(output_file, 'w') as outfile:

      for h,s in FastaReader(inputfile):
        stops = 9999
        translation = sixFrameTranslation(s)

        for frame in translation:
          st = translation[frame].count('*')
          if st < stops:
            best = frame
            stops = st

        if stops <= 2:
          outfile.write(">" + h + " frame_" + str(best) + "\n")
          outfile.write(translation[best] + "\n")
        else:
          badSeqs += 1
          bad_lengths.append(len(s))
          badfile.write(">" + h + "\n")
          badfile.write(s + "\n")

        seqCount += 1
  print ((100.0*badSeqs)/seqCount, "percent or ",
    badSeqs, " out of ", seqCount, " were not translated.")

  return output_file, output_file + "_BadSeqs"

def get_long_ORFS(translation, len_cutoff):
  keep = []
  for frame in translation:
    orfs = translation[frame].split('*')
    o_counts = 0
    for o in orfs:
      if len(o) >= len_cutoff:
        keep.append(("_frm"+str(frame)+"_orf"+str(o_counts), o))
  return keep

def pull_out_long_ORFs(bad_file, len_cutoff, outputdir):

  out_file = (outputdir + os.path.splitext(os.path.basename(bad_file))[0]
    + "_TranslongORFS.fa")

  with open(out_file, 'w') as outfile:
    for h,s in FastaReader(bad_file):
      translation = sixFrameTranslation(s)
      orfs = get_long_ORFS(translation, len_cutoff)
      for o in orfs:
        outfile.write(">" + h + o[0] + "\n")
        outfile.write(o[1] + "\n")

  return out_file

def searchhmmer(seqfile, hmmfile, bit_threshold, outdir, outname
    ,verbose , output_MA=False):
    #runs HMMERs search program
    path = HMMER3

    search_file = outdir +outname + "search.txt"

    if output_MA:
        ma_file = outdir + outname + "MA.sth"
        search_cmd = (path
            + " --nonull --noali --max"
            + " --incdomT " + str(bit_threshold)
            + " --domT " + str(bit_threshold)
            + " --nonull --noali --max"
            + " --incdomT " + str(bit_threshold)
            + " --domT " + str(bit_threshold)
            + " -A " + ma_file #multiple alignment output
            + " --domtblout " + search_file #hit output
            + " " + hmmfile
            + " " + seqfile
            + " > /dev/null" #throw away unwanted output
            )
    else:
        search_cmd = (path
            + " --nonull --noali --max"
            + " --incdomT " + str(bit_threshold)
            + " --domT " + str(bit_threshold)
            + " --domtblout " + search_file #hit output
            + " " + hmmfile
            + " " + seqfile
            + " > /dev/null") #throw away unwanted output
    if verbose:
        print search_cmd
    check_call(search_cmd, shell=True)

    if output_MA:
        return ma_file, search_file
    else:
        return search_file

def filter_with_HMMER(orfFile, hmmfile, outputdir):
  outname = os.path.splitext(os.path.basename(orfFile))[0]

  search_file = searchhmmer(orfFile, hmmfile, 0.1, outputdir, outname,True)

  output_file = (outputdir + os.path.splitext(os.path.basename(orfFile))[0]
    + "_matchedHMMER.fa")

  # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
  with open(search_file, 'r') as searchfile:
    for line in searchfile:
      if line[0]=='#':
        continue
      line = line.strip().split()
      target = line[0]
      query = line[3]

      print target, query

  return output_file


def merge_files(easy, hmmMatch, outputdir):

  out_file = (outputdir + os.path.splitext(os.path.basename(easy))[0]
    + "_finalTranslated.fa")

  with open(out_file, 'w') as outfile:
    for h,s in FastaReader(easy):
      outfile.write(">"+h+"\n")
      outfile.write(s+"\n")
    for h,s in FastaReader(hmmMatch):
      outfile.write(">"+h+"\n")
      outfile.write(s+"\n")

  print "Success!"

  return



def main():
  parser = OptionParser()

  parser.add_option("-c", "--contig", dest="contig_file",
      help="contig file of contigs to be translated")

  parser.add_option("", "--hmm", dest="hmmfile",
      help="hmmfile file of Rask block HMMs to be searched")

  parser.add_option("-l", "--lengthThreshold", dest="length", type=int
      , default=130
      , help="length of the shortest contig to be investigated")

  parser.add_option("-o", "--outputdir", dest="outputdir",
      help=("output directory. This will be created if it doesn't already"
          + " exist."))

  (options, args) = parser.parse_args()

  easy, hard = split_easy_from_hard(options.contig_file
    , options.outputdir)

  longORFs = pull_out_long_ORFs(hard, options.length, options.outputdir)

  hmmMatch = filter_with_HMMER(longORFs, options.hmmfile, options.outputdir)

  merge_files(easy, hmmMatch, options.outputdir)

  return

if __name__ == '__main__':
  main()

















