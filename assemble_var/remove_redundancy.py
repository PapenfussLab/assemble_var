import sys, os
from third_party_runners import run_blast
from mungo.fasta import FastaReader

OVERLAP = 0.7
IDENTITY = 95

sequences = sys.argv[1]
outdir = sys.argv[2]


#first run all-vs-all blast
bfile = run_blast(sequences, sequences, outdir, True)

#now get the lengths of all the contigs
contig_len = {}
for h,s in FastaReader(sequences):
    contig_len[h]=len(s)

#now iterate through search finding redundant contigs
redundant_contigs = set()
with open(bfile, 'r') as blastsearch:
    for line in blastsearch:
        tokens = line.strip().split()
        if tokens[0]==tokens[1]: #matching itself
            continue
        if float(tokens[2]) >= IDENTITY:
            # print "identity ",tokens[2]
            if ((float(tokens[3])/contig_len[tokens[0]] >= OVERLAP) or 
                (float(tokens[3])/contig_len[tokens[1]] >= OVERLAP)):
                # print "overlap", max(float(tokens[3])/contig_len[tokens[0]],float(tokens[3])/contig_len[tokens[1]])
                #we want to remove the shorted contig
                if contig_len[tokens[0]] < contig_len[tokens[1]]:
                    redundant_contigs.add(tokens[0])
                else:
                    redundant_contigs.add(tokens[1])

#now output a file without the redundant contigs
total_contigs = 0
count = 0
outname = os.path.splitext(os.path.basename(sequences))[0]
sample_name = "_".join(outname.split("_")[:2])
with open(outdir+outname+"_rmRedundant.fa", 'w') as outfile:
    for h,s in FastaReader(sequences):
        total_contigs+=1
        if h not in redundant_contigs:
            outfile.write(">" + sample_name + "_" + h + "\n")
            outfile.write(s+"\n")
        else:
            count+=1

print "inital contigs: ", total_contigs
print "number removed: ", count


