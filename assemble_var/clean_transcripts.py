from mungo.fasta import FastaReader
import sys

class Locus:
    def __init__(self, name):
        self.transcripts = []
        self.name = name
        self.max_trans_length = 0

class Transcript:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence



def remove_shorter_sequences(fasta_file, per_within_max, len_cutoff
    , outputfile):

    locus_dict = {}

    #Load in fasta file
    for h,s in FastaReader(fasta_file):
        locus_name = h.split("_")[1]
        if locus_name not in locus_dict:
            locus_dict[locus_name]=Locus(h)

        locus_dict[locus_name].transcripts.append(Transcript(h, s))

        if locus_dict[locus_name].max_trans_length < len(s):
            locus_dict[locus_name].max_trans_length = len(s)

    #now output with cleaning
    with open(outputfile,'w') as outfile:
        for locus_name in locus_dict:
            for t in locus_dict[locus_name].transcripts:
                if ((len(t.sequence)/float(locus_dict[locus_name].max_trans_length))
                    < (1-per_within_max)): #too far away from the longest
                    continue
                if len(t.sequence) < len_cutoff:
                    continue
                outfile.write(">"+t.name+"\n")
                outfile.write(t.sequence+"\n")





def main():
    remove_shorter_sequences(sys.argv[1], float(sys.argv[2])
        , float(sys.argv[3]), sys.argv[4])


if __name__ == '__main__':
    main()