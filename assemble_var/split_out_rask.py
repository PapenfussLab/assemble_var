import sys, os
from third_party_runners import run_blast
from mungo.fasta import FastaReader


contig_file = sys.argv[1]
blastdb_file = sys.argv[2]
outdir = sys.argv[3]


name = os.path.splitext(os.path.basename(contig_file))[0]

blast_out = run_blast(blastdb_file, contig_file, outdir, True)

rask_hits = set()
with open(blast_out , 'r') as blastfile:
    for line in blastfile:
        rask_hits.add(line.split()[0])

with open(outdir + name + "_rask.fa", 'w') as raskout:
    with open (outdir + name + "_nonrask.fa", 'w') as nonraskout:
        for h,s in FastaReader(contig_file):
            if h in rask_hits:
                raskout.write(">"+h+"\n")
                raskout.write(s+"\n")
            else:
                nonraskout.write(">"+h+"\n")
                nonraskout.write(s+"\n")


