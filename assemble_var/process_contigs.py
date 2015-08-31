import sys
import os
from optparse import OptionParser
from third_party_runners import run_blast, align_w_subread
from third_party_runners import convert_to_bam_create_index
from pysam import *
from mungo.fasta import FastaReader

HEADER = "hit id\tcontig id\tcontig length\tdata base name\tsequence id\talignment length\tpercent identity\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\trpk\n"

class Counter:
    def __init__(self):
        self.counts = 0
    def __call__(self, alignment):
        self.counts += 1

class BHit:
    def __init__(self, length, database_type, query_id, sub_id, per_identity, align_len, mismatches
        , gap_opens, qstart, qend, sstart, send, evalue, bitscore):
        self.query_id = query_id
        self.sub_id = sub_id
        self.per_identity = float(per_identity)
        self.align_len = int(align_len)
        self.mismatches = int(mismatches)
        self.gap_opens = int(gap_opens)
        self.query = [int(qstart),int(qend)]
        self.seq = [int(sstart), int(send)]
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.contig_length = int(length)
        self.database_type = database_type

        self.re_orientate()

        # print self.query, self.seq

    def re_orientate(self):
        if self.seq[0]>self.seq[1]:
            t = self.seq[1]
            self.seq[1] = self.seq[0]
            self.seq[0] = t
            #the contig is aligned in reverse
            self.query[0] = self.contig_length - self.query[1]
            self.query[1] = self.query[0]+self.align_len


        #error checking
        if self.seq[0]>=self.seq[1]:
            print self.query, self.seq
            raise NameError('sequence order problem (sequence)')
        if self.query[0]>=self.query[1]:
            print self.query, self.seq
            raise NameError('sequence order problem (query)')

        return

class Contig:
    def __init__(self, name, length, rpk=None):
        self.name = name
        self.length = int(length)
        self.rpk = rpk
        self.hits = []


def generate_summary(bamfile, outputdir, contig_file
    , blast_files=[]):
    blast_names = [os.path.splitext(os.path.basename(name))[0] for name in blast_files]

    #first create contig dictionary
    contig_dict = {}
    for h,s in FastaReader(contig_file):
        if len(h.split("_"))>1:
            contig_name = "_".join(h.split('_')[:4])
            print contig_name
        else:
            contig_name = h
        contig_dict[contig_name] = Contig(contig_name, len(s))

    #now add in the rpk info
    samfile = Samfile(bamfile, "rb")
    for h in samfile.references:
        c = Counter()
        samfile.fetch(h, callback=c)
        if len(h.split("_"))>1:
            print h
            h = "_".join(h.split('_')[:4])
            print h #rename to match dictionary
        else:
            h = h
        contig_dict[h].rpk = float(c.counts) / (contig_dict[h].length / 1e3)

    #Now add in the blast hit information
    for blast_file in blast_files:
        blast_name = os.path.splitext(os.path.basename(blast_file))[0]
        with open(blast_file, 'r') as bfile:
            for line in bfile:
                tokens = line.strip().split()
                if len(tokens[0].split("_"))>1:
                    name = "_".join(tokens[0].split('_')[:4])
                else:
                    name = tokens[0]
                contig_dict[name].hits.append(BHit(
                    contig_dict[name].length, blast_name
                    , *tokens))

    #Now write out the ouput
    # name_contig_file = os.path.splitext(os.path.basename(contig_file))[0]
    with open(outputdir + "analytics_output.txt", 'w') as outfile:
        outfile.write(HEADER)
        id_ = 0
        for contig in contig_dict:
            if not contig_dict[contig].hits: # it hasn't matched any of the databases
                outfile.write("\t".join([str(id_), contig_dict[contig].name, str(contig_dict[contig].length)]
                    + ["-"]*12 + [str(contig_dict[contig].rpk)])+"\n")
                id_+=1
            else:
                for hit in contig_dict[contig].hits:
                    outfile.write("\t".join([str(id_), contig_dict[contig].name, str(contig_dict[contig].length)]
                        + [hit.database_type, hit.sub_id]
                        + [str(hit.align_len), str(hit.per_identity)]
                        + [str(hit.mismatches), str(hit.gap_opens)]
                        + [str(hit.query[0]), str(hit.query[1])]
                        + [str(hit.seq[0]),str(hit.seq[1])]
                        + [str(hit.evalue), str(hit.bitscore)]
                        + [str(contig_dict[contig].rpk)]
                        )+"\n")
                    id_+=1

    return outputdir + "analytics_output.txt"


def analyse_contigs(contig_file, read1, read2, outputdir
   , fasta_ref_files=[], verbose=False):
    #first prepare blast files for analysis

    blast_files = []
    for reference in fasta_ref_files:
        blast_files.append(run_blast(reference, contig_file, outputdir
            , verbose))

    #now align reads to contigs
    samfile = align_w_subread(read1, read2, contig_file, outputdir, verbose)

    #now convert to bam
    bamfile = convert_to_bam_create_index(contig_file, samfile, verbose)

    # bamfile="notused"
    #now compute analytics
    outfile = generate_summary(bamfile, outputdir
        , contig_file, blast_files)

    return outfile

def trim_contigs(length, contig_file, outdir):
    out_file = outdir + "trim_transcripts.fa"
    with open(out_file, 'w') as outfile:
        for h,s in FastaReader(contig_file):
            if len(s)<length:
                continue
            outfile.write(">"+h+"\n")
            outfile.write(s+"\n")

    return out_file



def main():
    parser = OptionParser()
    parser.add_option("-r", "--read1", dest="read1",
        help="first set of read pairs")
    parser.add_option("-R", "--read2", dest="read2",
        help="second set of read pairs")

    parser.add_option("", '--ref',
                  type='string', action='append',
                  default=[], dest='ref_files',
                  help=('reference files that contain the sequence we want to '
                    + 'compare the contigs to'))
    parser.add_option("-o", "--outputdir", dest="outputdir",
        help="output directory")

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("-c", "--contig", dest="contig_file",
        help="contig file of contigs to be analysed")

    parser.add_option("", "--length", dest="length", type=int
        , help="length of the shortest contig to be investigated")


    (options, args) = parser.parse_args()

    #trim contigs that are less than length
    trim_contig = trim_contigs(options.length, options.contig_file
        , options.outputdir)

    #run if called from the command line
    analyse_contigs(trim_contig, options.read1, options.read2
        , options.outputdir, options.ref_files
        , options.verbose)


if __name__ == '__main__':
    main()