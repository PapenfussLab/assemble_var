import sys
from Bio import SeqIO
from numpy import mean
from numpy import median
from scipy.stats import *

"""
cut out transcripts that are shorter than min_transcript_len

only look at transcripts that are longer than the SMALLEST_LENGTH_FRACTION
compared to the longest transcript of the same locus

only keep loci with number of transcripts <= TperL
use "Inf" if no limit for TperL

print the excluded loci to the screen
"""

SMALLEST_LENGTH_FRACTION = 0.1

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python process_oases_transcripts.py oases_outputDIR max_TperL min_transcript_len outfile"
		sys.exit(0)
		
	DIR = sys.argv[1]+"/"
	if sys.argv[2] == "Inf":
		max_TperL = 999999
	else:
		max_TperL = int(sys.argv[2])
	min_transcript_len = int(sys.argv[3])
	out = sys.argv[4]
	
	#Get node length and coverage
	print ("reading stats")
	stats = open(DIR+"stats.txt","r")
	first = True #skip header line
	node_dict = {} #key is nodeid, value is average_coverage
	node_length = {} #key is nodeid, value is length
	for i in stats:
		if first == True:
			first = False
			continue
		spls = i.strip().split("\t")
		if spls[6] != 'Inf':
			node_dict[spls[0]] = float(spls[6])
	stats.close()
	
	print ("reading contig_ordering")
	c_ord = open(DIR+"contig-ordering.txt","r")
	count = 0
	locus_dict = {} #key is locus, value is list of transcript geomeans
	for i in SeqIO.parse(c_ord,"fasta"):
		if "Node" in i.id:
			node_id = i.id.split("_Node_")[1]
		if "Transcript" in i.id:
			locus_id = i.id.split("_Transcript_")[0]
			if locus_id not in locus_dict:
				locus_dict[locus_id] = []
			seqspls = str(i.seq).split("->")
			nodes = []
			for j in seqspls:
				tid = j.split(":")[0]
				if tid[0] == "-":
					tid = tid[1:]
				nodes.append(tid)
			nums = []
			transcript_size = 0
			for j in nodes:
				if j in node_dict: nums.append(node_dict[j])
			locus_dict[locus_id].append(gmean(nums))
		count += 1
	#remove loci with high numbers of transcripts
	removes = []
	for i in locus_dict:
		print i,locus_dict[i]
	for i in locus_dict:
		if len(locus_dict[i]) > max_TperL:
			removes.append(i)
	for i in removes:
		del locus_dict[i]
	
	locus_lengths = {}#key is locus, value is list of transcript lengths
	locus_longest = {}#key is locus, value is the length of longest transcript
	trans = open(DIR+"transcripts.fa","r")
	for i in SeqIO.parse(trans,"fasta"):
		locus_id = i.id.split("_Transcript_")[0]
		transcript_size = len(i.seq)
		if locus_id not in locus_lengths:
			locus_lengths[locus_id] = []
		locus_lengths[locus_id].append(transcript_size)
		if locus_id not in locus_longest:
			locus_longest[locus_id] = 0
		if locus_longest[locus_id] < transcript_size:
			locus_longest[locus_id] = transcript_size
	

	#do the smallest_length_fraction here
	best_trans_num = {}#key is locus, value is the index of the best transcript
					   #(n+1 is what will be in the transcripts file)
	best_trans_val = {}#key is locus, value is the value of the coverage
	for i in locus_dict:
		tlist = locus_dict[i]#coverage
		llist = locus_lengths[i]#length
		besttr = 0
		number_of_best = 0
		for number in range(len(tlist)):
			if tlist[number] > besttr and ((llist[number]/float(locus_longest[i]))>SMALLEST_LENGTH_FRACTION):
				besttr = tlist[number]
				number_of_best = number
		best_trans_num[i] = number_of_best
		best_trans_val[i] = besttr
	
	print ("writing scatterplot")
	ScatterPlot = open(DIR+'scatterplot_mine.csv',"w")
	ScatterPlotString='Locus\tTranscript\tPercentLen\tCoverageFrac'
	ScatterPlot.write(ScatterPlotString+'\n')
	for i in locus_dict:
		tlist = locus_dict[i]
		llist = locus_lengths[i]
		for number in range(len(tlist)):
			PercentLen='%.2f' % ((llist[number]*100/float(locus_longest[i])))
			if tlist[number]/best_trans_val[i] > 1.0:
				# max(tlist) > 2.0:
				CoverageFrac=0.0 #ones that didn't pass the length
			else:
				CoverageFrac='%.4f' % float(tlist[number]/best_trans_val[i])
			ScatterPlotString=str(str(i.split("Locus_")[1])+'\t'+str(number+1)+'\t'+PercentLen+'\t'+str(CoverageFrac))
			ScatterPlot.write(ScatterPlotString+'\n')
	
	print ("getting best transcript")
	#do the absolute length cutoff here
	keep_seqs = []
	trans.seek(0,0)
	for i in SeqIO.parse(trans,"fasta"):
		if "Transcript" in i.id:
			locus_id = i.id.split("_Transcript_")[0]
			try:
				tran_id = best_trans_num[locus_id]
				if i.id.split("_Transcript_")[1].split("/")[0]  == str(tran_id+1):
					if len(i.seq) > min_transcript_len:
						keep_seqs.append(i)
			except:#too many transcripts
				continue
	trans.close()
	
	print ("writing transcripts to : "+out)
	outfile = open(out,"w")
	SeqIO.write(keep_seqs,outfile,"fasta")
	outfile.close()
	
