#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>


#include "subread.h"
#include "input-files.h"
#include "core.h"

static struct option long_options[] =
{
	{"index",  required_argument, 0, 'i'},
	{"read",  required_argument, 0, 'r'},
	{"read2",  required_argument, 0, 'R'},
	{"output",  required_argument, 0, 'o'},
	{"subreads",  required_argument, 0, 'n'},
	{"threads",  required_argument, 0, 'T'},
	{"indel",  required_argument, 0, 'I'},
	{"phred",  required_argument, 0, 'P'},
	{"minmatch",  required_argument, 0, 'm'},
	{"minmatch2",  required_argument, 0, 'p'},
	{"mindist",  required_argument, 0, 'd'},
	{"maxdist",  required_argument, 0, 'D'},
	{"order",  required_argument, 0, 'S'},
	{"DPMismatch",  required_argument, 0, 'X'},
	{"DPMatch",  required_argument, 0, 'Y'},
	{"DPGapOpen",  required_argument, 0, 'G'},
	{"DPGapExt",  required_argument, 0, 'E'},
	{"junction", no_argument, 0, 'J'},
	{"unique",  no_argument, 0, 'u'},
	{"color-convert",  no_argument, 0, 'b'},
	{"multi",  required_argument, 0, 'B'},
	{"hamming",  no_argument, 0, 'H'},
	{"quality",  no_argument, 0, 'Q'},
	{"trim5", required_argument, 0, '5'},
	{"trim3", required_argument, 0, '3'},
	{"memoryMultiplex",  required_argument, 0, 0},
	{"rg",  required_argument, 0, 0},
	{"rg-id",  required_argument, 0, 0},
	{"BAMoutput", no_argument, 0, 0},
	{"BAMinput", no_argument, 0, 0},
	{"SAMinput", no_argument, 0, 0},
	{"reportPairedMultiBest",  no_argument, 0, 0},
	{"reportFusions", no_argument, 0, 0},
	{"gzFASTQinput", no_argument, 0, 0},
	{"extraColumns",  no_argument, 0, 0},
	{"forcedPE",  no_argument, 0, 0},
	{"ignoreUnmapped",  no_argument, 0, 0},
	{"accurateFusions",  no_argument, 0, 0},
	{"maxMismatches",  required_argument, 0, 'M'},
	{0, 0, 0, 0}
};


void print_usage_core_aligner()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);
	SUBREADputs("Usage:");
	SUBREADputs("");
	SUBREADputs(" ./subread-align [options] -i <index_name> -r <input> -o <output>");
	SUBREADputs("");
	SUBREADputs("Required arguments:");
	SUBREADputs("    ");
	SUBREADputs("    -i --index     <index>  base name of the index.");
	SUBREADputs("   ");
	SUBREADputs("    -r --read      <input>  name of the input file(FASTQ/FASTA format by default");
	SUBREADputs("                            . See below for more supported formats). Both base-");
	SUBREADputs("                            space and color-space read data are supported. For");
	SUBREADputs("                            paired-end reads, this gives the first read file");
	SUBREADputs("                            and the other read file should be specified using");
	SUBREADputs("                            the -R option.");
	SUBREADputs("    ");
	SUBREADputs("Optional general arguments:");
	SUBREADputs("    ");
	SUBREADputs("    -o --output    <output> name of the output file(SAM format by default). If");
	SUBREADputs("                            not provided, mapping results will be output to the");
	SUBREADputs("                            standard output (stdout).");
	SUBREADputs("");
	SUBREADputs("    -n --subreads  <int>    number of selected subreads, 10 by default.");
	SUBREADputs("    ");
	SUBREADputs("    -m --minmatch  <int>    consensus threshold (minimal number of consensus");
	SUBREADputs("                            subreads required) for reporting a hit. If paired-");
	SUBREADputs("                            end read data are provided, this gives the consensus");
	SUBREADputs("                            threshold for the read which receives more votes");
	SUBREADputs("                            than the other read from the same pair. 3 by default");
	SUBREADputs("    ");
	SUBREADputs("    -T --threads   <int>    number of threads, 1 by default.");
	SUBREADputs("    ");
	SUBREADputs("    -I --indel     <int>    number of indels allowed, 5 by default. Indels of up");
	SUBREADputs("                            to 200bp long can be detected.");
	SUBREADputs("    ");
	SUBREADputs("    -B --multi     <int>    Specify the maximal number of equally-best mapping");
	SUBREADputs("                            locations allowed to be reported for each read. 1");
	SUBREADputs("                            by default. Allowed values are between 1 to 16");
	SUBREADputs("                            (inclusive). 'NH' tag is used to indicate how many");
	SUBREADputs("                            alignments are reported for the read and 'HI' tag");
	SUBREADputs("                            is used for numbering the alignments reported for");
	SUBREADputs("                            the same read, in the output. Note that -u option");
	SUBREADputs("                            takes precedence over -B.");
	SUBREADputs("");
	SUBREADputs("    -P --phred     <3:6>    the format of Phred scores in input files, '3' for");
	SUBREADputs("                            phred+33 and '6' for phred+64. '3' by default.");
	SUBREADputs("");
	SUBREADputs("    -u --unique             only uniquely mapped reads will be reported (reads");
	SUBREADputs("                            mapped to multiple locations in the reference genome");
	SUBREADputs("                            will not be reported). This option can be used");
	SUBREADputs("                            together with option '-H' or '-Q'.");
	SUBREADputs("");
	SUBREADputs("    -Q --quality            using mapping quality scores to break ties when more");
	SUBREADputs("                            than one best mapping location is found.");
	SUBREADputs("");
	SUBREADputs("    -H --hamming            using Hamming distance to break ties when more than");
	SUBREADputs("                            one best mapping location is found.");
	SUBREADputs("");
	SUBREADputs("    -b --color-convert      convert color-space read bases to base-space read");
	SUBREADputs("                            bases in the mapping output. Note that the mapping");
	SUBREADputs("                            itself will still be performed at color-space.");
	SUBREADputs("");
	SUBREADputs("    -M --maxMismatches <int> Specify the maximum number of mis-matched bases");
	SUBREADputs("                            allowed in the alignment. 3 by default. Mis-matches");
	SUBREADputs("                            found in soft-clipped bases are not counted.");
	SUBREADputs("   ");
	SUBREADputs("       --reportFusions      report discovered genomic fusion events such as");
	SUBREADputs("                            chimeras. Discovered fusions will be saved to a file");
	SUBREADputs("                            (*.fusions.txt). Detailed mapping results for fusion");
	SUBREADputs("                            reads will be saved to the SAM/BAM output file as");
	SUBREADputs("                            well. Secondary alignments of fusion reads will be");
	SUBREADputs("                            saved to the following optional fields: CC(Chr),");
	SUBREADputs("                            CP(Position), CG(CIGAR) and CT(strand). Note that");
	SUBREADputs("                            each fusion read occupies only one row in the");
	SUBREADputs("                            SAM/BAM output file.");
	SUBREADputs("");
	SUBREADputs("       --trim5     <int>    trim off <int> number of bases from 5' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --trim3     <int>    trim off <int> number of bases from 3' end of each");
	SUBREADputs("                            read. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --rg-id     <string> specify the read group ID. If specified,the read");
	SUBREADputs("                            group ID will be added to the read group header");
	SUBREADputs("                            field and also to each read in the mapping output.");
	SUBREADputs("");
	SUBREADputs("       --rg        <string> add a <tag:value> to the read group (RG) header in");
	SUBREADputs("                            in the mapping output.");
	SUBREADputs("");
	SUBREADputs("       --gzFASTQinput       specify that the input read data is in gzipped");
	SUBREADputs("                            FASTQ/FASTA format.");
	SUBREADputs("");
	SUBREADputs("       --SAMinput           specify that the input read data is in SAM format.");
	SUBREADputs("");
	SUBREADputs("       --BAMinput           specify that the input read data is in BAM format.");
	SUBREADputs("");
	SUBREADputs("       --BAMoutput          specify that mapping results are saved into a BAM");
	SUBREADputs("                            format file.");
	SUBREADputs("");
	SUBREADputs("       --DPGapOpen  <int>   a numeric value giving the penalty for opening a");
	SUBREADputs("                            gap when using the Smith-Waterman dynamic");
	SUBREADputs("                            programming algorithm to detect insertions and");
	SUBREADputs("                            deletions. The Smith-Waterman algorithm is only");
	SUBREADputs("                            applied for those reads which are found to contain");
	SUBREADputs("                            insertions or deletions. -1 by default.");
	SUBREADputs("");
	SUBREADputs("       --DPGapExt   <int>   a numeric value giving the penalty for extending the");
	SUBREADputs("                            gap, used by the Smith-Waterman algorithm. 0 by");
	SUBREADputs("                            default.");
	SUBREADputs("");
	SUBREADputs("       --DPMismatch <int>   a numeric value giving the penalty for mismatches,");
	SUBREADputs("                            used by the Smith-Waterman algorithm. 0 by default.");
	SUBREADputs("");
	SUBREADputs("       --DPMatch    <int>   a numeric value giving the score for matches used by");
	SUBREADputs("                            the Smith-Waterman algorithm. 2 by default.");
	SUBREADputs("");
	SUBREADputs("    -v                      output version of the program.");
	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs("Optional arguments for paired-end reads:");
	SUBREADputs("");
	SUBREADputs("    -R --read2     <input>  name of the second input file. The program will then");
	SUBREADputs("                            be switched to the paired-end read mapping mode.");
	SUBREADputs("");
	SUBREADputs("    -p --minmatch2 <int>    consensus threshold for the read which receives less");
	SUBREADputs("                            votes than the other read from the same pair, 1 by");
	SUBREADputs("                            default.");
	SUBREADputs("");
	SUBREADputs("    -d --mindist   <int>    minimum fragment/template length, 50bp by default.");
	SUBREADputs("");
	SUBREADputs("    -D --maxdist   <int>    maximum fragment/template length, 600bp by default.");
	SUBREADputs("");
	SUBREADputs("    -S --order     <ff:fr:rf> orientation of the two read from the same pair,");
	SUBREADputs("                            'fr' by default.");
	SUBREADputs("");
	SUBREADputs("");
	SUBREADputs("For more information about these arguments, please refer to the User Manual.");
	SUBREADputs("");

}


int parse_opts_aligner(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	
	int is_64_bit_computer = sizeof(char *)>4; 

	optind = 0;
	opterr = 1;
	optopt = 63;

	global_context->config.entry_program_name = CORE_PROGRAM_SUBREAD;
	global_context->config.max_mismatch_exonic_reads = 3;
	global_context->config.max_mismatch_junction_reads = 3;
	global_context->config.use_dynamic_programming_indel = 1;

	// config.extending_search_indels is changed from 1 to 0 on 10/mar/2014
	global_context->config.extending_search_indels = 0;
	global_context->config.big_margin_record_size = 9; 

	if(argc<2)
	{
		print_usage_core_aligner();
		return -1;
	}

/*

	for(c = 0; c<argc; c++)
	{
		printf("[%d]\t\t%s\n", c , argv[c]);
	}

*/

	while ((c = getopt_long (argc, argv, "xsvJS:L:AHd:D:n:m:p:G:E:X:Y:P:R:r:i:l:o:T:I:t:B:bFcuUfM:Q1:2:3:5:?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'v':
				core_version_number("Subread-align");
				return -1;
			case 'G':
				global_context->config.DP_penalty_create_gap = atoi(optarg);
				break;
			case 'Y':
				global_context->config.DP_match_score = atoi(optarg);
				break;
			case 'E':
				global_context->config.DP_penalty_extend_gap = atoi(optarg);
				break;
			case 'X':
				global_context->config.DP_mismatch_penalty = atoi(optarg);
				break;
			case '3':
				global_context->config.read_trim_3 = atoi(optarg); 
				break;
			case '5':
				global_context->config.read_trim_5 = atoi(optarg); 
				break;
			case 'J':
				global_context->config.show_soft_cliping = 0;
				break;
			case 'B':
				global_context->config.multi_best_reads = atoi(optarg); 
				if(global_context->config.multi_best_reads<1)
					global_context->config.multi_best_reads=1;
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[1]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 'b':
				global_context->config.convert_color_to_base = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				//global_context->config.total_subreads = min(31,global_context->config.total_subreads );
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atof(optarg);
				break;
			case 'T':
				global_context->config.all_threads = atoi(optarg);
				if(global_context->config.all_threads <1) global_context->config.all_threads = 1;
				if(global_context->config.all_threads >32) global_context->config.all_threads = 32;

				break;
			case 'r':
				strncpy(global_context->config.first_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'R':
				global_context->input_reads.is_paired_end_reads = 1;
				strncpy(global_context->config.second_read_file, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'i':
				strncpy(global_context->config.index_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'o':
				strncpy(global_context->config.output_prefix, optarg, MAX_FILE_NAME_LENGTH-1);
				break;
			case 'I':
				global_context->config.max_indel_length = atoi(optarg);

				if(!is_64_bit_computer) global_context->config.max_indel_length = min(global_context->config.max_indel_length , 16); 
				if(global_context->config.max_indel_length <0)global_context->config.max_indel_length =0;
				if(global_context->config.max_indel_length > MAX_INSERTION_LENGTH)global_context->config.max_indel_length = MAX_INSERTION_LENGTH;
				if(global_context->config.max_indel_length > 16)
				{
					global_context->config.reassembly_subread_length = 16;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 3;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 3;
					global_context->config.reassembly_key_length = 28;

					global_context->config.is_third_iteration_running = 1;


					// These options were removed for maximising the sensitivity.
					// They should be put back when we want a higher indel accuracy.
					//global_context->config.max_mismatch_exonic_reads = 1;
					//global_context->config.max_mismatch_junction_reads = 1;
					//global_context->config.total_subreads = 28;
					//global_context->config.minimum_subread_for_first_read = 3;
					//global_context->config.minimum_subread_for_second_read = 1;
					//global_context->config.do_big_margin_filtering_for_reads = 0;
					//global_context->config.extending_search_indels = 0;
					//global_context->config.use_dynamic_programming_indel = 1;
					//global_context->config.flanking_subread_indel_mismatch = 0;

					global_context->config.do_superlong_indel_detection = 0;
				}
				break;
			case 'M':
				global_context->config.max_mismatch_exonic_reads = atoi(optarg);
				break;
			case 'P':
				if (optarg[0]=='3')
					global_context->config.phred_score_format = FASTQ_PHRED33;
				else
					global_context->config.phred_score_format = FASTQ_PHRED64;
				break;
			case 'p':
				global_context->config.minimum_subread_for_second_read = atoi(optarg);
				break;
			case 't':
				sprintf(global_context->config.temp_file_prefix, "%s/core-temp-sum-%06u-%05u", optarg, getpid(), rand());
				break;
			case 'F':
				global_context->config.is_second_iteration_running = 0;
				global_context->config.report_sam_file = 0;
				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
			case 'Q':
				global_context->config.use_quality_score_break_ties = 1;
				break;
				
			case 0:
				if(strcmp("memoryMultiplex", long_options[option_index].name)==0) 
				{
					global_context->config.memory_use_multiplex = atof(optarg);
				}
				else if(strcmp("rg-id", long_options[option_index].name)==0) 
				{
					strcpy(global_context->config.read_group_id, optarg);
				}
				else if(strcmp("rg", long_options[option_index].name)==0) 
				{
					strcat(global_context->config.read_group_txt, "\t");
					strcat(global_context->config.read_group_txt, optarg);
				}
				else if(strcmp("BAMoutput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_output = 1;
				}
				else if(strcmp("BAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 1;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("gzFASTQinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_gzip_fastq=1;
				}
				else if(strcmp("extraColumns", long_options[option_index].name)==0) 
				{
					global_context->config.SAM_extra_columns=1;
				}
				else if(strcmp("SAMinput", long_options[option_index].name)==0) 
				{
					global_context->config.is_BAM_input = 0;
					global_context->config.is_SAM_file_input = 1;
				}
				else if(strcmp("reportPairedMultiBest", long_options[option_index].name)==0) 
				{
					global_context->config.report_multiple_best_in_pairs = 1;
				}
				else if(strcmp("ignoreUnmapped", long_options[option_index].name)==0) 
				{
					global_context->config.ignore_unmapped_reads = 1;
				}
				else if(strcmp("reportFusions", long_options[option_index].name)==0) 
				{
					global_context->config.is_rna_seq_reads = 1;
					global_context->config.do_fusion_detection = 1;
					global_context->config.prefer_donor_receptor_junctions = 0;
					global_context->config.do_big_margin_filtering_for_reads = 1;
				}
				else if(strcmp("accurateFusions", long_options[option_index].name)==0) 
				{
					global_context->config.more_accurate_fusions = 1;
				}
				break;
			case '?':
			default:
				SUBREADprintf("Unknown option: -%c",c);
				print_usage_core_aligner();
				return -1 ;
		}
	}

	global_context->config.more_accurate_fusions = global_context->config.more_accurate_fusions && global_context->config.do_fusion_detection;
	if(global_context->config.more_accurate_fusions)
	{
		global_context->config.high_quality_base_threshold = 999999;
		global_context->config.max_mismatch_junction_reads = 0;
		global_context->config.do_big_margin_filtering_for_junctions = 1;
		global_context->config.total_subreads = 20;
	}

	if(global_context->config.is_SAM_file_input) global_context->config.phred_score_format = FASTQ_PHRED33;

	return 0;
}





#if defined MAKE_STANDALONE
int main(int argc , char ** argv)
{
#elif defined RUNNING_ENV_JAVA
int subread_aligner_main(int argc , char ** argv)
{
#else
int main_align(int argc , char ** argv)
{
#endif

//	printf("SIZE_OF_ALN=%d\n", sizeof(alignment_result_t));
//	printf("SIZE_OF_VOT=%d\n", sizeof(voting_context_t));
	return core_main(argc, argv, parse_opts_aligner);
}

