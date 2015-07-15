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
  
  
#ifndef SUBREAD_CORE_INDEL_H_
#define SUBREAD_CORE_INDEL_H_

#include "subread.h"
#include "hashtable.h"
#include "core.h"

// chromosome events can be indels, junctions or fusions.
// if it is an insertion event, event_large_site = event_small_site+1.

//#define MAX_EVENT_ENTRIES_PER_SITE 5
#define MAX_EVENT_ENTRIES_PER_SITE 12 
#define CHRO_EVENT_TYPE_REMOVED 0
#define CHRO_EVENT_TYPE_INDEL 8
#define CHRO_EVENT_TYPE_LONG_INDEL 16 
#define CHRO_EVENT_TYPE_POTENTIAL_INDEL 32 
#define CHRO_EVENT_TYPE_JUNCTION 64 
#define CHRO_EVENT_TYPE_FUSION 128
#define CHRO_EVENT_TYPE_SNP 256 

#define EVENT_SEARCH_BY_SMALL_SIDE 10
#define EVENT_SEARCH_BY_LARGE_SIDE 20
#define EVENT_SEARCH_BY_BOTH_SIDES 30


#define REASSEMBLY_WINDOW_LENGTH 350

//#define is_target_window_X(x) ((x + 1) * REASSEMBLY_WINDOW_LENGTH / 2 >= (10734463 % BASE_BLOCK_LENGTH) && (x- 1) * REASSEMBLY_WINDOW_LENGTH /2-1 <= (10734463%BASE_BLOCK_LENGTH) )
#define is_target_window_X(x) 0
//#define MAXIMUM_EVENT_NUMBER 300000

typedef struct{
	unsigned int event_small_side;
	unsigned int event_large_side;
	//union
	//{
		short indel_length;
		short junction_flanking_left;
	//};
	short junction_flanking_right;

	unsigned char event_type;
	char indel_at_junction;
	char is_negative_strand;	// this only works to junction detection, according to 'GT/AG' or 'CT/AC' donors. This only applys to junctions.
	char is_strand_jumped;		// "strand jumped" means that the left and right sides are on different strands. This only applys to fusions.
	char is_donor_found;		// only for junctions: GT/AG is found at the location.
						// Also, if "is_strand_jumped" is true, all coordinates (e.g., splicing points, cover_start, cover_end, etc) are on "reversed read" view.

	//char is_ambiguous;
	char connected_next_event_distance;	// the distance (negative or positive) to the next event in the table. For example, if the cigar string is 10M3I1M1I10M, event "3I" will have 1 here .
	char connected_previous_event_distance;	// the distance (negative or positive) to the next event in the table. For example, if the cigar string is 10M3I1M1I10M, event "1I" will have 1 here.

	//char inserted_bases[(1+MAX_INSERTION_LENGTH) / 4 + 1];
	char * inserted_bases;
	unsigned short supporting_reads;
	unsigned short anti_supporting_reads;
	unsigned short final_counted_reads;
	unsigned short final_reads_mismatches;
	unsigned int global_event_id;
	float event_quality;
} chromosome_event_t;


struct reassmebly_window_allele
{
	char rebuilt_window[8000];
	float allele_quality;
	int rebuilt_size;
};

typedef struct{
	gehash_t * voting_indexes;
	char * chro_name;
	unsigned long long int * start_keys;
	short * start_offsets;

	unsigned int * read_no_counter;
	unsigned int block_start_linear_pos;
	HashTable * read_sequence_table;
	HashTable * read_position_table;
	HashTable * read_quality_table;
	gene_vote_t * vote_list;
	gene_vote_t * vote_list_rectify;
	short * read_rectify_space;

	char rebuilt_window[2500];
	int rebuilt_window_size;


	struct reassmebly_window_allele * final_alleles;

	unsigned int used_read_ids[2000];
	int used_read_number;


	int search_cost;
	int total_matched_bases;
	int max_matched_bases;
	unsigned int window_start_pos;
} reassembly_by_voting_block_context_t;



typedef struct{
	HashTable ** de_bruijn_graphs;
	char * chro_name;
	unsigned long long int * start_keys;
	short * start_offsets;

	unsigned int block_start_linear_pos;
} reassembly_block_context_t;


typedef struct{
	HashTable * event_entry_table;
	unsigned int total_events;
	unsigned int current_max_event_number;
	chromosome_event_t * event_space_dynamic;
	HashTable * local_reassembly_pileup_files;

	short ** dynamic_align_table;
	char ** dynamic_align_table_mask;
} indel_context_t;

typedef struct{
	HashTable * event_entry_table;
	unsigned int total_events;
	unsigned int current_max_event_number;
	chromosome_event_t * event_space_dynamic;
	unsigned short * final_counted_reads_array;
	unsigned short * final_reads_mismatches_array;

	short ** dynamic_align_table;
	char ** dynamic_align_table_mask;
} indel_thread_context_t;

int init_indel_tables(global_context_t * context);
int destroy_indel_module(global_context_t * context);
int init_indel_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task);
int finalise_indel_thread(global_context_t * global_context, thread_context_t * thread_context, int task);
int find_new_indels(global_context_t * global_context, thread_context_t * thread_context, int pair_number, char * read_name, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id);
int write_indel_final_results(global_context_t * context);
int search_event(global_context_t * global_context,HashTable * event_table, chromosome_event_t * event_space, unsigned int pos, int search_type, char event_type, chromosome_event_t ** return_buffer);

void set_alignment_result(global_context_t * global_context, int pair_number, int is_second_read, int best_read_id, unsigned int position, int votes, gene_vote_number_t * indel_record, short best_cover_start, short best_cover_end, int is_negative_strand, unsigned int minor_position, unsigned int minor_votes, unsigned int minor_coverage_start, unsigned int minor_coverage_end, unsigned int split_point, int inserted_bases, int is_strand_jumped, int is_GT_AG_donors, int used_subreads_in_vote, int noninformative_subreads_in_vote, int major_indel_offset, int minor_indel_offset);

void put_new_event(HashTable * event_table, chromosome_event_t * new_event , int event_no);
void remove_neighbour(global_context_t * global_context);
int build_local_reassembly(global_context_t *global_context , thread_context_t *thread_context , int pair_number, char * read_name_1 , char * read_text_1 ,char * qual_text_1 , int read_len_1, int read_len_2, int is_second_read, int best_read_id, int is_paired_unmapped);
int finalise_long_insertions(global_context_t * global_context);

// This function sets the global context with default values.
void init_global_context(global_context_t * context);

int write_local_reassembly(global_context_t *global_context, HashTable *pileup_fp_table, unsigned int anchor_pos, char * read_name , char * read_text ,char * qual_text , int read_len, int is_anchor_certain);

int finalise_long_insertions_by_hashtable(global_context_t * global_context);

void destroy_pileup_table(HashTable* local_reassembly_pileup_files);

chromosome_event_t * reallocate_event_space(global_context_t* global_context,thread_context_t* thread_context,int event_no);

int there_are_events_in_range(char * bitmap, unsigned int pos, int sec_len);

int anti_supporting_read_scan(global_context_t * global_context);
#endif
