import sys, os
from mungo.fasta import FastaReader
from collections import defaultdict
from optparse import OptionParser

from third_party_runners import run_blast


def get_longest(contig_file, outfile, min_length):
    #returns a dictionary of header sequence pairs of the longest contig in 
    #each locus
    locus_dict = defaultdict(list)
    contig_dict = {}
    for h,s in FastaReader(contig_file):
        if len(h.split("_"))>1:
            locus_dict[h.split("_")[1]].append((h,s))
        else:
            contig_dict[h] = s

    longest_contigs = []

    for locus in locus_dict:
        max_len = 0
        for h,s in locus_dict[locus]:
            if len(s)>max_len:
                max_len = len(s)
                curr = (h,s)
        longest_contigs.append(curr)

    with open(outfile, 'w') as outfas:
        for contig in longest_contigs:
            if len(contig[1]) >= min_length:
                outfas.write(">" + contig[0] + "\n")
                outfas.write(contig[1] + "\n")
        for contig in contig_dict:
            if len(contig_dict[contig]) >= min_length:
                outfas.write(">" + contig + "\n")
                outfas.write(contig_dict[contig] + "\n")

    return longest_contigs


def filter_ref_with_blast(fasta_ref_files, contig_file, percent_overlap
    , outfile, outdir):

    #first get list of contigs
    contigs = {}
    for h,s in FastaReader(contig_file):
        contigs[h] = s

    print "Number of contigs before filtering: ", len(contigs.keys())

    #now run blast against the reference files which we want not to be 
    #present in the data i.e. human
    blast_files = []
    for reference in fasta_ref_files:
        blast_files.append(run_blast(reference, contig_file, outdir
            , True))

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
    for name in bad_contigs:
        del contigs[name]

    #now write resulting contigs to a file
    with open(outfile, 'w') as outfas:
        for contig in contigs:
            outfas.write(">" + contig + "\n")
            outfas.write(contigs[contig] + "\n")

    print "Number of contigs after filtering: ", len(contigs.keys())


def filter_contigs(contig_file, fasta_ref_files, percent_overlap
    , outdir, min_length):
    name = os.path.splitext(os.path.basename(contig_file))[0]
    get_longest(contig_file, outdir + name + "_longest.fas", min_length)

    filter_ref_with_blast(fasta_ref_files, outdir + name + "_longest.fas"
        , percent_overlap, outdir + name + "_filtered.fas", outdir)



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


    (options, args) = parser.parse_args()

    filter_contigs(options.contig_file, options.ref_files, options.percent
        , options.outputdir, options.min_length)


if __name__ == '__main__':
    main()

