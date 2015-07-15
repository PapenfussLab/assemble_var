from optparse import OptionParser
import third_party_runners as thrd
import os


def build(options):
    #create the output directory if it doesn't exist
    try:
        os.mkdir(options.outputdir)
    except OSError, e:
        if e.errno != 17: #ignores error if folders has already been created.
            raise
        pass


    #first create a temporary directory in the output direectory
    try:
        outputdir = options.outputdir + "temp_files/"
        os.mkdir(outputdir)
    except OSError, e:
        if e.errno != 17: #ignores error if folders has already been created.
            raise
        pass

    curr_dir = os.getcwd()
    os.chdir(outputdir)
    
    # # print thrd.getScriptPath()
    # read1, read2 = thrd.trim_galore(options.read1, options.read2, outputdir
    #     , options.verbose)

    # alignment = thrd.align_w_subread(read1, read2, options.reference, outputdir
    #     , options.verbose, options.ref_index)

    # bamfile = thrd.convert_to_bam_create_index(options.reference, alignment
    #     , options.verbose)

    # single_reads, paired_reads = thrd.extract_possible_var_reads(bamfile
    #     , outputdir, options.var_files, options.verbose)

    # # single_reads, paired_reads = outputdir + "extracted_unpaired.fq", outputdir + "extracted_paired.fq"

    # if options.pear:
    #     single_reads, paired_reads = thrd.merge_with_pear(single_reads
    #         , paired_reads, outputdir, options.verbose)

    # #now digital normalisation on the single_reads and paired reads seperately
    # if options.norm:
    #     single_reads, paired_reads = thrd.digi_norm(single_reads
    #         , paired_reads, outputdir, options.verbose)



    single_reads = "/home/users/allstaff/tonkin-hill.g/novel_data/assemble_all_renormalised/temp_files/re_normalised_single.fa"
    paired_reads = "/home/users/allstaff/tonkin-hill.g/novel_data/assemble_all_renormalised/temp_files/re_normalised_paired.fa"

    # transcript_file = thrd.assemble_paired_reads(single_reads, paired_reads
    #     , options.outputdir, options.ins_length, options.verbose)
    
    folder_prefix = thrd.assemble_paired_reads(single_reads, paired_reads
        , options.outputdir, options.ins_length, options.verbose)
    
    transcript_file_61, harsh_transcipt_file = thrd.filter_Locus_1(options.outputdir
        , folder_prefix, options.verbose)

    thrd.assemble_contigs_cap3(harsh_transcipt_file, options.outputdir
        , verbose=options.verbose
        , outname="harsh_transcipt_file_cap3.fa")

    thrd.assemble_contigs_cap3(transcript_file_61, options.outputdir
        , verbose=options.verbose
        , outname="transcript_file_61_cap3.fa")

    os.chdir(curr_dir)



def main():
    parser = OptionParser()
    parser.add_option("-r", "--read1", dest="read1",
        help="first set of read pairs")
    parser.add_option("-R", "--read2", dest="read2",
        help="second set of read pairs")
    parser.add_option("", "--reference", dest="reference"
        , help=("a fasta file containing the reference"
            + " genomes we want to filter out"))
    parser.add_option("", "--index", dest="ref_index", default=False,
        help=("the location of the index files."
            + " Usefule if performing multiple assemblies as this"
            + " does not need to be recomputed each time."))

    parser.add_option('-v', '--var',
                  type='string', action='append',
                  default=[], dest='var_files',
                  help=('var files that contain var regions we want to keep'
                    + ' - that is a subset of reference'))

    parser.add_option("-o", "--outputdir", dest="outputdir",
        help="output directory")

    parser.add_option("", "--ins_length", dest="ins_length"
        , default=False, help="the insert length passed to oases")

    parser.add_option("","--verbose", action="store_true", dest="verbose"
        , default=False, help="turns on more detailed output")

    parser.add_option("", "--pear", action="store_true", dest="pear"
        , default=False, help="merge read pairs that overlap before oases.")

    parser.add_option("", "--norm", action="store_true", dest="norm"
        , default=False, help="merge read pairs that overlap before oases.")

    # #optional parts of the assembly are called like this
    # parser.add_option(,"--nofilter", action="store_false", dest=isfilter, default=True
    #     , help="turns off filtering of non-var plasmodium falciparum sequence")

    (options, args) = parser.parse_args()

    build(options)

if __name__ == '__main__':
    main()