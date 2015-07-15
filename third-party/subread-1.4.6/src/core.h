/***************************************************************

   The Subread software package is free software package: 
   you can redistribute it and/or modify it under the terms
   of the GNU General Public License as published by the 
   Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Subread is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   See the GNU General Public License for more details.

   Authors: Drs Yang Liao and Wei Shi

  ***************************************************************/
  
  
//Core.c is the totally refactored body of the current subread, subjunc, subfusion and subindel programs.
//This program runs in this way:
// 1, it loads reads from FASTQ, FASTA, PLAIN or SAM files
// 2, it change the read and the quality strings according to FF/FR/RR/RF orders, such that every read is in the "F" manner.
//    Namely, the second read in a pair is reversed.
// 3, it builds voting table for each read.
// 4, it performs the first iteration.
// 5, it performs the second iteration.
// 6, it reports results regarding the reads. 
// 7, if the input files contain too many reads (e.g., greater than 14M in total), step 1 ~ 6 are repeated until all reads are finalised.
// 8, it reports global results.

//Core.c maintains the following context:
// 1, the global context.
// 2, the overall context for each module.
// 3, the local context for voting each read for each module.
// The structure of each context is named as
//     Context_ModuleName_Scope

//Thread context is maintained for each thread.

// There should be a function for each module to store/load its global context.

#ifndef _SUBREAD_CORE_H_
#define _SUBREAD_CORE_H_
#include <pthread.h> 
#include "subread.h"
#include "sambam-file.h"

//#define is_target(s) (memcmp((s),"chr901_242237_242737_0:0:0_0:0:0_3761", 22)==0)

#define is_target(s) 0

//#define _global_retrieve_voting_context(global_context, pair_number) (global_context->chunk_vote_records[pair_number])

#define _global_retrieve_alignment(global_context, pair_number, is_second_read, best_read_id) ((global_context->input_reads.is_paired_end_reads?(global_context -> chunk_alignment_records[(2*pair_number+is_second_read)* global_context->config.multi_best_reads + best_read_id]):(global_context -> chunk_alignment_records[pair_number*global_context->config.multi_best_reads + best_read_id])))

#define _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id) ((global_context->input_reads.is_paired_end_reads?(global_context -> chunk_alignment_records + (2*pair_number+is_second_read)* global_context->config.multi_best_reads + best_read_id):(global_context -> chunk_alignment_records+pair_number*global_context->config.multi_best_reads + best_read_id)))

#define _global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, best_read_id) ((global_context->input_reads.is_paired_end_reads?(global_context -> chunk_subjunc_records + (2*pair_number+is_second_read)* global_context->config.multi_best_reads + best_read_id):(global_context -> chunk_subjunc_records+pair_number*global_context->config.multi_best_reads + best_read_id)))


#define _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read)  (is_second_read?(global_context->big_margin_record + (2*pair_number+is_second_read)* global_context->config.big_margin_record_size):(global_context->big_margin_record + pair_number * global_context->config.big_margin_record_size))

#define mark_gapped_read(res) (res)-> result_flags|= CORE_IS_GAPPED_READ;

typedef struct{
	int is_paired_end_reads;
	gene_input_t first_read_file;
	gene_input_t second_read_file;
	unsigned long long first_read_file_size;
	unsigned long long first_file_blocks[64];
	unsigned long long second_file_blocks[64];
	unsigned int reads_in_blocks[64];
	unsigned int start_read_number_blocks[64];
	double avg_read_length;
} read_input_t;



typedef struct{
	// running_scheme
	int all_threads;
	int is_first_iteration_running;
	int is_second_iteration_running;
	int is_third_iteration_running;
	float memory_use_multiplex;
	char temp_file_prefix[MAX_FILE_NAME_LENGTH];
	unsigned int reads_per_chunk;

	// input_scheme
	char first_read_file[MAX_FILE_NAME_LENGTH];
	char second_read_file[MAX_FILE_NAME_LENGTH];
	char medium_result_prefix[MAX_FILE_NAME_LENGTH];


	short read_trim_5;
	short read_trim_3;
	int is_first_read_reversed;
	int is_second_read_reversed;
	int space_type;
	int is_methylation_reads;
	int phred_score_format;
	int is_SAM_file_input;
	int is_gzip_fastq;


	// reporting scheme
	char read_group_id[MAX_FILE_NAME_LENGTH];
	char read_group_txt[MAX_FILE_NAME_LENGTH];
	char output_prefix[MAX_FILE_NAME_LENGTH];
	int report_sam_file;
	int report_no_unpaired_reads;
	int min_mapped_fraction;
	int max_mismatch_exonic_reads;
	int max_mismatch_junction_reads;
	int ignore_unmapped_reads;
	int report_unmapped_using_mate_pos;
	int report_multi_mapping_reads;
	int downscale_mapping_quality;
	int is_BAM_input;
	int is_BAM_output;
	int convert_color_to_base;
	int SAM_extra_columns;
	int report_multiple_best_in_pairs;
	unsigned int multi_best_reads;

	// basic voting
	char index_prefix[MAX_FILE_NAME_LENGTH];
	int total_subreads;
	int minimum_subread_for_first_read;
	int minimum_subread_for_second_read;
	float minimum_exonic_subread_fraction;
	int minimum_pair_distance;
	int maximum_pair_distance;
	int restrected_read_order;
	int show_soft_cliping;
	int max_indel_length;
	int expected_pair_distance;
	int ambiguous_mapping_tolerance;
	int use_hamming_distance_break_ties;
	int use_quality_score_break_ties;
	int big_margin_record_size;

	// subjunc
	int entry_program_name;
	char is_rna_seq_reads;
	char do_big_margin_filtering_for_junctions;
	char do_big_margin_filtering_for_reads;
	char limited_tree_scan;
	char use_hamming_distance_in_exon;
	unsigned int maximum_intron_length;
	int high_quality_base_threshold;
	char max_insertion_at_junctions;
	char check_donor_at_junctions;

	// subfusion
	int do_fusion_detection;
	int prefer_donor_receptor_junctions;
	int more_accurate_fusions;

	// indel
	char do_superlong_indel_detection;
	char extending_search_indels;
	int k_mer_length;
	int reassembly_start_read_number;
	int reassembly_key_length;
	int reassembly_subread_length;
	int reassembly_window_multiplex;
	int reassembly_tolerable_voting;
	int reassembly_window_alleles;
	int init_max_event_number;
	int use_dynamic_programming_indel;
	int use_bitmap_event_table;
	int flanking_subread_indel_mismatch;
	int DP_penalty_create_gap;
	int DP_penalty_extend_gap;
	int DP_match_score;
	int DP_mismatch_penalty;

} configuration_t;

#define CORE_IS_GT_AG_DONORS 1
#define CORE_NOTFOUND_DONORS 2
#define CORE_IS_STRAND_JUMPED 4
#define CORE_IS_NEGATIVE_STRAND 8
#define CORE_IS_FULLY_EXPLAINED 16
#define CORE_IS_BREAKEVEN 32
#define CORE_IS_GAPPED_READ 64 

#define CORE_CIGAR_OPT_M 0
#define CORE_CIGAR_OPT_S 1
#define CORE_CIGAR_OPT_D 2
#define CORE_CIGAR_OPT_I 3
#define CORE_CIGAR_OPT_B 4
#define CORE_CIGAR_OPT_N 5
#define CORE_CIGAR_OPT_BB 6
#define CORE_CIGAR_OPT_NN 7

#define CORE_PROGRAM_SUBREAD 100
#define CORE_PROGRAM_SUBJUNC 200
#define CORE_PROGRAM_SUBINDEL 1000

typedef struct
{

	unsigned int selected_position;
	// 4 bytes
	gene_vote_number_t selected_votes;
	gene_vote_number_t used_subreads_in_vote;
	unsigned char noninformative_subreads_in_vote;
	unsigned char final_quality;
	// this coverage is the range on reads, in point of view of "main piece" strand (i.e., "is_negative_strand")
	char indels_in_confident_coverage; 
	char result_flags;
	// 12 bytes
	union{
		struct{
			gene_vote_number_t selected_indel_record [MAX_INDEL_SECTIONS*3 + 1];
			unsigned short confident_coverage_start; 
			unsigned short confident_coverage_end; 
		};
		char cigar_string[MAX_INDEL_SECTIONS * 3+5];
	};
	// 4x bytes

	union
	{
		unsigned long long Score_L;
		struct{
			short final_mismatched_bases;
			short best_second_diff_bases;
		};
	};
	// 48 bytes
	unsigned long long int Score_H;
	// 56 butes

} alignment_result_t;

#define CORE_MAX_CIGAR_LEN (MAX_INDEL_SECTIONS * 3+5)
#define CORE_MAX_CIGAR_STR_LEN 110
#define CORE_ADDITIONAL_INFO_LENGTH 400

typedef struct
{
	short split_point;
	gene_vote_number_t minor_votes;
	char  double_indel_offset;
	char indel_at_junction;
	unsigned int minor_position;
	// this coverage is the range on reads, in point of view of "main piece" strand
	unsigned short minor_coverage_start; 
	unsigned short minor_coverage_end; 

} subjunc_result_t;

typedef struct{
	int thread_id;
	pthread_t thread;

	// module_thread_context is different to module_context, though some contents in module_thread_context could be also in module context (e.g., junction tables)
	// modules functions may deside to use objects in which context.
	void * module_thread_contexts[5];
	gene_value_index_t * current_value_index;
	unsigned int processed_reads_in_chunk;

	// per chunk parameters
	gene_input_t * ginp1;
	gene_input_t * ginp2;
	unsigned int reads_to_be_done;
	unsigned int read_block_start;
} thread_context_t;


typedef struct{
	// basic running configuration
	configuration_t config;

	// for the index
	gehash_t * current_index;
	gene_value_index_t * current_value_index;
	gene_value_index_t all_value_indexes[100];
	int index_block_number;
	int current_index_block_number;
	int will_remove_input_file;
	int is_phred_warning;

	// global locks
	subread_lock_t thread_initial_lock;

	// global FILE*
	SamBam_Writer * output_bam_writer;
	FILE * output_sam_fp;
	FILE * long_insertion_FASTA_fp;

	// running contexts
	void * module_contexts[5]; 
	read_input_t input_reads;
	alignment_result_t * chunk_alignment_records;	// arrangement: PE:: array_offset = ( read_pair_no * 2 + is_second_read ) * best_read_number + best_read_id
							// arrangement: SE:: array_offset = read_pair_no * best_read_number + best_read_id
	subjunc_result_t * chunk_subjunc_records;	// arrangement: PE:: array_offset = ( read_pair_no * 2 + is_second_read ) * best_read_number + best_read_id
							// arrangement: SE:: array_offset = read_pair_no * best_read_number + best_read_id
	unsigned char * big_margin_record;
	gene_offset_t chromosome_table;
	double start_time;
	double align_start_time;

	unsigned long long all_processed_reads;
	unsigned long long all_mapped_reads;
	unsigned long long all_correct_PE_reads;
	unsigned int all_junctions;
	unsigned int all_fusions;
	unsigned int all_indels;

	unsigned long long current_circle_start_position_file1;
	unsigned long long current_circle_start_position_file2;
	unsigned long long current_circle_end_position_file1;
	unsigned long long current_circle_end_position_file2;
	unsigned int processed_reads_in_chunk;

	// per chunk parameters
	unsigned int reads_to_be_done;
	unsigned int read_block_start;

} global_context_t;


#define MODULE_INDEL_ID 0
#define MODULE_JUNCTION_ID 1

#define STEP_VOTING 10
#define STEP_ITERATION_ONE 20
#define STEP_ITERATION_TWO 30
#define STEP_ITERATION_THREE 40
#define STEP_WRITE_CHUNK_RESULTS 90

#define MEMORY_OPTIMISATION_LAPTOP 10
#define MEMORY_OPTIMISATION_DESKTOP 20
#define MEMORY_OPTIMISATION_SUPERCOMPUTER 100

// This function initialise the module contexts that are needed according to the running configuration.
int init_modules(global_context_t * context);

// This function undo any memory allocation and close opened files in init_modules.
int destroy_modules(global_context_t * context);

// This function shows welcome messages and shows the current configuration to stdout.
// This function also checks the legality of the parameters in context->config
int print_configuration(global_context_t * context);

// This function opens all related files except the index. It however counts the piece number of the index.
int load_global_context(global_context_t * context);

// This function deallocates memory and closes files that were allowcated and opened in load_global_context.
int destroy_global_context(global_context_t * context);


// This function processes the read chunks (24M reads each), including three steps: 1, voting; 2, iteration one; 3, iteration two.
int read_chunk_circles(global_context_t * context);


// write the final results in each module. Note that the SAM file is writen in read_chunk_circles, after each read chunk is processed.
int write_final_results(global_context_t * context);

// this function is called by an interface to start the main process with a customized opt parser.
int core_main(int argc , char ** argv, int (parse_opts (int , char **, global_context_t * )));

// decompress a cigar string to cigar
// cigar_len is the maximum length of decompressed cigar
// bincigar is the maximum length of compressed cigar
// it returns the length of decompressed cigar, or -1 if buffer is too short.
int bincigar2cigar(char * cigar, int cigar_len, char * bincigar, int bincigar_max_len, int read_len);

// compress a cigar string to bincigar
// bincigar_len is the maximum length of compressed cigar
// it returns the length of compressed cigar, or -1 if buffer is too short.
int cigar2bincigar(char *cigar, char *bincigar, int bincigar_len);

// print the logo of our program
void print_subread_logo();

// print a line in the box
void print_in_box(int line_width, int is_boundary, int is_center, char * pattern,...);

// find the value index covering this read
// it returns NULL if no index is found.
gene_value_index_t * find_current_value_index(global_context_t * global_context, unsigned int pos, int len);

// generate the time string 
void char_strftime(char * tbuf);

int term_strncpy(char * dst, char * src, int max_dst_memory);

int is_result_in_PE(alignment_result_t * aln);

void core_version_number(char * program);


// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar);

int chimeric_cigar_parts(global_context_t * global_context , unsigned int sel_pos, int is_first_section_negative_strand, int is_first_section_reversed, char * in_cigar, unsigned int * out_poses, char ** out_cigars, char * out_strands, int read_len, short * out_read_lens);

void warning_file_limit();
void quick_sort(void * arr,int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

// L_Minus_R should return -1, 0 or 1 when L<R, L==R or L>R.
// The result is from Small to Large.
void merge_sort(void * arr, int arr_size, int L_Minus_R (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge_SmallFirst(void * arr, int start, int items, int items2));
#endif
