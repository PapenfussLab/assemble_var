import os, sys, inspect
from subprocess import check_call

# add path to mungo library
cmd_subfolder = (os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    + "/third-party/")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import re
import glob
from mungo.fasta import FastaReader
from collections import defaultdict

def getScriptPath():
    directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    return directory

def trim_galore(inputfile, inputfile2, outputdir, adapter, verbose=False):
    #filter out ILLUMINA primer etc
    #returns the filenames of the filtered reads
    scriptPath = getScriptPath()

    trim_cmd = (scriptPath
        + "/third-party/trim_galore/trim_galore "
        + " --phred33 "
        + " --paired "
        + " --dont_gzip")

    # if adapter:
    #     trim_cmd = trim_cmd + " --adapter " + adapter

    trim_cmd = (trim_cmd
        + " --output_dir " + outputdir + " "
        + inputfile + " "
        + inputfile2
        )

    if verbose:
        print trim_cmd
    # check_call(trim_cmd, shell=True)


    return (outputdir + os.path.splitext(os.path.basename(inputfile))[0]+"_val_1.fq"
        , outputdir + os.path.splitext(os.path.basename(inputfile2))[0]+"_val_2.fq")

#TODO: split out indexing command so its only run once and clean up where its sending the indexes!!
def align_w_subread(read1, read2, reference, outputdir, verbose=False
    , index=False, phred=False):
    #aligns reads quickly to the reference genome
    #returns alignment filename
    scriptPath = getScriptPath()
    scriptPath = scriptPath +"/third-party/subread-1.4.6/bin/"

    if not index:
        index_cmd = (scriptPath
            +"subread-buildindex "
            + "-M 12000 "
            + "-o ref_index "
            + reference)
        index = "ref_index"
        if verbose:
            print index_cmd
        #first build indices for subread program
        # check_call(index_cmd, shell=True)



    align_cmd = (scriptPath
        + "subread-align "
        + " -d 50"
        + " -D 600"
        + " -i " + index
        + " -r " + read1
        + " -R " + read2
        + " -o " + outputdir + "alignment.sam")

    if phred:
        align_cmd = align_cmd + " --phred " + str(phred)

    if verbose:
        print align_cmd
    #now align read to reference
    # check_call(align_cmd, shell=True)

    return outputdir + "alignment.sam"

def convert_to_bam_create_index(reference, filename, verbose=False):
    #converts filename into a bam with the reference and puts the result into outputdir

    sam_cmd = ("samtools view -b"
        + " -T " + reference
        + " -o " + ".".join(filename.split(".")[:-1]) + ".bam"
        + " " + filename)
    sort_cmd = ("samtools sort "
        + ".".join(filename.split(".")[:-1]) + ".bam "
        + ".".join(filename.split(".")[:-1]) + "_sort")
    index_cmd = ("samtools index "
        + ".".join(filename.split(".")[:-1]) + "_sort.bam "
        # + ".".join(filename.split(".")[:-1]) + "_sort"
        )

    if verbose:
        print sam_cmd
        print sort_cmd
        print index_cmd

    # check_call(sam_cmd, shell=True)
    # check_call(sort_cmd, shell=True)
    # check_call(index_cmd, shell=True)

    return ".".join(filename.split(".")[:-1]) + "_sort.bam"

def extract_possible_var_reads(bamfile, outputdir, varfiles=None, verbose=False):
    #first we take out the reads that have not been aligned.
    scriptPath = getScriptPath()

    count = 1
    count_var = 1

    # for f, F in zip([4,8,12], [264,260,256]): #used if we want mate pairs as well
    for f in [4]: #just get unmapped reads
        sam_cmd = ("samtools view -u"
            + " -f " + str(f)
            # + " -F " + str(F)
            + " " + bamfile
            + " | samtools sort -n - " + outputdir + "temp_file"+str(count)
            # + " > " + outputdir + "temp_file"+str(count)+".bam"
            )
        if verbose:
            print sam_cmd
        # check_call(sam_cmd, shell=True)
        count+=1

    #now extract reads that have aligned to the var regions
    for varfile in varfiles:
        bed_cmd = (scriptPath + "/third-party/bedtools2/bin/bedtools intersect -ubam"
            + " -a " + bamfile
            + " -b " + varfile
            + " | samtools sort -n - " + outputdir + "temp_file"+str(count)
            # + " > " + outputdir + "temp_file"+str(count_var)+".sam"
            )
        if verbose:
            print bed_cmd
        # check_call(bed_cmd, shell=True)
        count+=1

    if len(varfiles)<2:
        varname = "temp_file1.sam"
    else:
        varname = "temp_file["+"".join([str(i) for i in range(1,count_var)])+"].sam"

    #now merge files together ready for assembly
    merge_cmd = ("samtools merge -f -u - "
        + outputdir + "temp_file["+"".join([str(i) for i in range(1,count)])+"].bam "
        # + outputdir + varname
        + " | samtools sort -n - " + outputdir + "extracted"
        )
    sort_cmd = ("samtools sort -n "
        + outputdir + "extracted.bam "
        + outputdir + "extracted_n")

    bamtofastq_cmd = (scriptPath + "/third-party/bamUtil_1.0.12/bamUtil/bin/bam bam2fastq"
        + " --merge"
        + " --in " + outputdir + "extracted_n" +".bam"
        + " --firstOut " + outputdir + "extracted_paired.fq"
        + " --unpairedOut " + outputdir + "extracted_unpaired.fq")

    if verbose:
        print merge_cmd
        print sort_cmd
        print bamtofastq_cmd
    # check_call(merge_cmd, shell=True)
    # check_call(sort_cmd, shell=True)
    # check_call(bamtofastq_cmd, shell = True)

    return outputdir + "extracted_unpaired.fq", outputdir + "extracted_paired.fq"

def merge_with_pear(single_seq, paired_seq, outputdir, verbose=False):
    scriptPath = getScriptPath()
    script = "python " + scriptPath + "/third-party/khmer/scripts/split-paired-reads.py"
    #first split files so PEAR can use them
    split_cmd = (script
        + " " + paired_seq)
    if verbose:
        print split_cmd

    # check_call(split_cmd, shell=True)

    #now run PEAR
    outname = ".".join(paired_seq.split(".")[:-1]) + "PEAR"
    pear_cmd = (scriptPath + "/third-party/pear"
        + " -f " + paired_seq + ".1"
        + " -r " + paired_seq + ".2"
        + " -o " + outname)
    if verbose:
        print pear_cmd

    # check_call(pear_cmd, shell=True)

    #now merge paired reads into single file
    script = "python " + scriptPath + "/third-party/khmer/scripts/interleave-reads.py"
    merge_cmd = (script
        + " " + outname + ".unassembled.forward.fastq"
        + " " + outname + ".unassembled.reverse.fastq"
        + " -o " + outputdir + "pairedPear.fastq"
        )
    if verbose:
        print merge_cmd
    # check_call(merge_cmd, shell=True)

    #Now merge the signle reads together
    cat_cmd = ("cat " + single_seq
        + " " + outname + ".assembled.fastq"
        + " > " + outputdir + "singlePear.fastq")
    if verbose:
        print cat_cmd
    # check_call(cat_cmd, shell=True)

    return outputdir + "singlePear.fastq", outputdir + "pairedPear.fastq"


# def split_out_pe(fastafile, verbose=False):
#     scriptPath = getScriptPath()
#     scriptPath = scriptPath + "/third-party/khmer/scripts/"
#     split_cmd = ("python "
#         + scriptPath + "extract-paired-reads.py "
#         + fastafile)

#     if verbose:
#         print split_cmd
#     check_call(split_cmd, shell=True)

#     return fastafile+".se", fastafile+".pe"

def assemble_paired_reads(fasta_single, fasta_paired, outputdir
    , ins_length=False, verbose=False):
    # print "WARNING: using default insert length which is not reliable according to oases manual"


    #de-novo assemble reads using oases
    scriptPath = getScriptPath()
    velveth = scriptPath + "/third-party/velveth"
    velvetg = scriptPath + "/third-party/velvetg"
    oases = scriptPath + "/third-party/oases/oases"
    scriptPath = scriptPath + "/third-party/oases/scripts/"


    assemble_cmd = ("python "
        + scriptPath + "oases_pipeline.py"
        + " -m 21 -M 65 -s 4"
        + " --velveth " + velveth
        + " --velvetg " + velvetg
        + " --oases " + oases
        + " --single" #only run single k assemblies
        + " -o assembly"
        + " -d \" -fastq -short " + fasta_single
        + " -shortPaired " + fasta_paired + " \"")
    if ins_length:
        assemble_cmd = assemble_cmd + " -p \" -ins_length " + ins_length + " \""


    if verbose:
        print assemble_cmd
    check_call(assemble_cmd, shell=True)
    check_call("mv assembly* " + outputdir, shell=True)


    # return outputdir + "assemblyMerged/transcripts.fa"
    return "assembly"

def filter_Locus_1(outputdir, folder_prefix, verbose=False):
    final_transcripts = []
    locus_dict = defaultdict(list)
    for transcript_file in glob.glob(outputdir + folder_prefix + '*/transcripts.fa'):
        #now need to extract all transcripts that are the only member of
        #their Locus
        name = transcript_file.strip('/transcripts.fa')[-2:]
        print name
        print transcript_file
        for h,s in FastaReader(transcript_file):
            locus = h.split("_")[1] + "_" + name
            locus_dict[locus].append((h + "_K" + name, s))
    #Now we want to write out all Locus' of length 1 to the outfile
    out_fileA = outputdir + "filter_locus1_transcripts_keep61.fa"
    with open(out_fileA, 'w') as outfile:
        for locus in locus_dict:
            if (locus.split("_")[-1]=='61') or (len(locus_dict[locus]) == 1):
                for t in locus_dict[locus]:
                    outfile.write(">" + t[0] + "\n")
                    outfile.write(t[1] + "\n")

    # out_fileB = outputdir + "filter_locus1_transcripts.fa"
    # with open(out_fileB, 'w') as outfile:
    #     for locus in locus_dict:
    #         if len(locus_dict[locus]) == 1:
    #             for t in locus_dict[locus]:
    #                 outfile.write(">" + t[0] + "\n")
    #                 outfile.write(t[1] + "\n")

    return out_fileA

def assemble_contigs_cap3(transcript_file, outputdir, verbose=False
    , outname="transcripts_cap3.fa"):

    #assemble the contigs using cap3
    cap3_cmd = ("cap3 "
        + transcript_file
        + " -p 99 -c 200"
        + " > " + outputdir + "cap.out")

    cat_cmd = ("cat "
        + transcript_file + ".cap.singlets "
        + transcript_file + ".cap.contigs "
        + "> " + outputdir + outname)

    if verbose:
        print cap3_cmd
        print cat_cmd
    check_call(cap3_cmd, shell=True)
    check_call(cat_cmd, shell=True)

    return outputdir + outname

def run_blast(reference, query_contigs, outdir, verbose=False):
    #first build blast data base
    #check if blastDB exists
    blastDB = (outdir + "blastDB" + "_"
        + os.path.splitext(os.path.basename(reference))[0])
    if os.path.isfile(blastDB + ".nhr") or os.path.isfile(blastDB + ".00.nhr"):
        if verbose:
            print "Blast DB already present continuing.."
    else:
        db_cmd = ("makeblastdb "
            + "-dbtype nucl "
            + "-in " + reference + " "
            + " -out " + blastDB)

        if verbose:
            print db_cmd
        check_call(db_cmd, shell=True)

    #now query contigs against DB
    output = (outdir + "blastHits" + "_"
        + os.path.splitext(os.path.basename(reference))[0] + "_"
        + os.path.splitext(os.path.basename(query_contigs))[0]
        )
    blast_cmd = ("blastn "
        + "-evalue 1e-5 -outfmt 6 "
        # + "-num_alignments " + str(num_hits) + " "
        + "-db " + blastDB + " "
        + "-query " + query_contigs + " "
        + "-out " + output
        )
    if verbose:
        print blast_cmd

    check_call(blast_cmd, shell=True)

    return output

def  digi_norm(single, paired, outputdir, verbose=False):
    scriptPath = getScriptPath()
    script = "python " + scriptPath + "/third-party/khmer/scripts/normalize-by-median.py"

    norm_cmd = (script
        + " -C 20 -k 20 -N 4 -x 2e9"
        + " -o " + outputdir + "normalised_single.fa"
        + " " + single)
    if verbose:
        print norm_cmd
    # check_call(norm_cmd, shell=True)

    norm_cmd = (script
        + " -C 20 -k 20 -N 4 -x 2e9"
        + " -p"
        + " -o " + outputdir + "normalised_paired.fa"
        + " " + paired)
    if verbose:
        print norm_cmd
    # check_call(norm_cmd, shell=True)

    return outputdir + "normalised_single.fq", outputdir + "normalised_paired.fq"


def combine_paired(read1, read2, outputdir, verbose=False):
    #function if we skip the alignment step.
    scriptPath = getScriptPath()
    script = "python " + scriptPath + "/third-party/khmer/scripts/interleave-reads.py"

    combine_cmd = (script
        + " -o " + outputdir + "interleaved_paired.fq"
        + " " + read1
        + " " + read2)

    if verbose:
        print combine_cmd

    check_call(combine_cmd, shell=True)

    open(outputdir + "interleaved_single.fq", 'a').close()

    return outputdir + "interleaved_single.fq", outputdir + "interleaved_paired.fq"


def assemble_paired_reads_soapDeNovoTrans(fasta_single, fasta_paired, outputdir
    , ins_length=False, verbose=False):

    scriptPath = getScriptPath()
    script = "python " + scriptPath + "/third-party/khmer/scripts/split-paired-reads.py"
    #first split files so PEAR can use them
    split_cmd = (script
        + " " + fasta_paired)

    if verbose:
        print split_cmd
    check_call(split_cmd, shell=True)

    config_str =(
"""#maximal read length
max_rd_len=250
[LIB]
#maximal read length in this lib
rd_len_cutof=250
#average insert size
avg_ins=0
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fasta file for single reads
q1=""" + fasta_paired + ".1\n"
+ "q2=" + fasta_paired + ".2\n"
+ "q=" + fasta_single)


    scriptPath = getScriptPath()
    script = (scriptPath
        + "/third-party/SOAPdenovo-Trans-bin-v1.03/SOAPdenovo-Trans-127mer")

    for K in ["21", "31", "41", "51", "61"]:
        config_file = outputdir + "soapConfigK" + K
        outputGraph = outputdir + "soapGraphK" + K

        soap_cmd = (script +
            " all"
            + " -s " + config_file
            + " -o " + outputGraph
            + " -K " + K
            + " -L 50" #not actually important as we only use the contigs
            )

        with open(config_file, 'w') as outfile:
            outfile.write(config_str)

        if verbose:
            print soap_cmd

        check_call(soap_cmd, shell=True)

    #now to merge with cap3
    # check_call("cat ")

    # assemble_contigs_cap3(transcript_file, outputdir, verbose=False
    #     , outname="transcripts_cap3.fa"):

    return








