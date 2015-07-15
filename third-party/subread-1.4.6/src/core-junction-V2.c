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
  
  
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "subread.h"
#include "sublog.h"
#include "gene-value-index.h"
#include "gene-algorithms.h"
#include "input-files.h"
#include "core.h"
#include "core-indel.h"
#include "core-junction.h"

int localPointerCmp_forbed(const void *pointer1, const void *pointer2)
{
	paired_exon_key *p1 = (paired_exon_key *)pointer1;
	paired_exon_key *p2 = (paired_exon_key *)pointer2;
	return !((p1-> big_key == p2 -> big_key) && (p2-> small_key == p1-> small_key));
}

unsigned long localPointerHashFunction_forbed(const void *pointer)
{
	paired_exon_key *p  = (paired_exon_key *)pointer;
	return p-> big_key ^ p-> small_key  ^ (p->big_key>> 15);
}

int localPointerCmp_forpos(const void *pointer1, const void *pointer2)
{
	return pointer1 != pointer2;
}

unsigned long localPointerHashFunction_forpos(const void *pointer)
{

	return (unsigned long) pointer & 0xffffffff;
}


typedef struct{
	unsigned int piece_main_abs_offset;
	unsigned int piece_minor_abs_offset;
	int piece_main_masks;
	short piece_main_coverage_start;
	short piece_main_coverage_end;

	short piece_main_hamming_match;
	short piece_main_read_quality;
	short piece_minor_hamming_match;
	short piece_minor_read_quality;
	short intron_length;

	char  *piece_main_indel_record;
	short piece_main_indels;
	short piece_minor_indel_offset;
	unsigned char piece_main_votes;
	unsigned char piece_minor_votes;

	short piece_minor_coverage_start;
	short piece_minor_coverage_end;
	short split_point;
	char is_GT_AG_donors;
	char is_donor_found;
	char is_strand_jumped;

	unsigned long long int Score_H;
	unsigned int Score_L;
} select_junction_record_t;


// read_head_abs_pos is the offset of the FIRST WANTED base.
void search_events_to_front(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_head_abs_offset, short remainder_len, short sofar_matched)
{
	short tested_read_pos;

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}


	int event_search_method;
	if(global_context -> config.do_fusion_detection)
		event_search_method = EVENT_SEARCH_BY_BOTH_SIDES;
	else
		event_search_method = EVENT_SEARCH_BY_SMALL_SIDE;

	// tested_read_pos is the index of the first base unwanted!
	if(MAX_EVENTS_IN_READ - 1> explain_context -> tmp_search_sections)
		for(tested_read_pos = 1 + 15 ; tested_read_pos <= remainder_len; tested_read_pos++)
		{
			int xk1, matched_bases_to_site;
			chromosome_event_t *site_events[MAX_EVENT_ENTRIES_PER_SITE+1];

			int jump_penalty = 0;

			unsigned potential_event_pos;
			if(explain_context -> current_is_strand_jumped)
				potential_event_pos = read_head_abs_offset - tested_read_pos +1;
			else
				potential_event_pos = read_head_abs_offset + tested_read_pos -1;
			int site_events_no = search_event(global_context, event_table , event_space , potential_event_pos, event_search_method , CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION , site_events);
			/*if(memcmp(explain_context->read_name, "HKJMKOB02G7RDV",14) == 0)
			{
				printf("FOUND THE EVENT FRONT:%d at %u\n", site_events_no, potential_event_pos);
				if(site_events_no)
					printf("EVENT0_type = %d\n", site_events[0]->event_type);
			}*/

			//if(explain_context -> pair_number==2074) printf("FF OFFSET=%d; LEDGE=%u; FOUND=%d\n", tested_read_pos, potential_event_pos, site_events_no);
			if(!site_events_no)continue;

			unsigned int tested_chro_begin;
			if(explain_context -> current_is_strand_jumped)
				tested_chro_begin = read_head_abs_offset - tested_read_pos + 1;
			else
				tested_chro_begin = read_head_abs_offset;

			matched_bases_to_site = match_chro(read_text, value_index, tested_chro_begin , tested_read_pos, explain_context -> current_is_strand_jumped, global_context -> config.space_type);


			//if(memcmp(explain_context->read_name, "HKJMKOB02G7RDV",14) == 0)
			//	printf("JUMP?%d > %d\n", (1+matched_bases_to_site)*10000 / tested_read_pos , 9000);

			if((1+matched_bases_to_site)*10000/tested_read_pos > 9000)
				for(xk1 = 0; xk1 < site_events_no ; xk1++)
				{
					chromosome_event_t * tested_event = site_events[xk1];

					// note that these two values are the index of the first wanted base.
					unsigned int new_read_head_abs_offset;

					if(global_context -> config.do_fusion_detection && tested_event->event_type == CHRO_EVENT_TYPE_FUSION)
						new_read_head_abs_offset = (potential_event_pos == tested_event -> event_large_side)?tested_event -> event_small_side:tested_event -> event_large_side;
					else
						new_read_head_abs_offset = tested_event -> event_large_side;


					short new_remainder_len = remainder_len - tested_read_pos + min(0, tested_event->indel_length);

					//int is_ambiguous = tested_event -> is_ambiguous;

					if(new_remainder_len>0)// && (new_remainder_len>8 || !is_ambiguous))
					{
						//if(explain_context -> pair_number==2074) printf("JUMPPED IN!\n");
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_end = explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start + tested_read_pos;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].event_after_section = tested_event;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].is_connected_to_large_side = (potential_event_pos == tested_event -> event_large_side);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].read_pos_start = tested_read_pos - min(0, tested_event -> indel_length);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].abs_offset_for_start = new_read_head_abs_offset;
						explain_context -> tmp_jump_length += (tested_event->event_large_side - tested_event->event_small_side);

						//if(tested_event->event_type == CHRO_EVENT_TYPE_FUSION) jump_penalty = 1;
						//else if(tested_event->event_type == CHRO_EVENT_TYPE_JUNCTION) jump_penalty = 1;

						int current_is_jumped = explain_context -> current_is_strand_jumped ;
						if(tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> is_strand_jumped)
							explain_context -> current_is_strand_jumped = !explain_context -> current_is_strand_jumped;

						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_strand_jumped = explain_context -> current_is_strand_jumped;

						explain_context -> tmp_search_sections ++;
						search_events_to_front(global_context, thread_context, explain_context, read_text + tested_read_pos -  min(0, tested_event->indel_length), qual_text + tested_read_pos -  min(0, tested_event->indel_length), new_read_head_abs_offset, new_remainder_len, sofar_matched + matched_bases_to_site - jump_penalty);
						explain_context -> tmp_search_sections --;

						explain_context -> current_is_strand_jumped = current_is_jumped;
						explain_context -> tmp_jump_length -= (tested_event->event_large_side - tested_event->event_small_side);
					}
					//if(global_context ->config.limited_tree_scan) break;
				}
			if(global_context ->config.limited_tree_scan && explain_context -> full_read_len <= EXON_LONG_READ_LENGTH) break;
		}

	int whole_section_matched = match_chro(read_text , value_index, explain_context -> current_is_strand_jumped?read_head_abs_offset - remainder_len +1:read_head_abs_offset, remainder_len , explain_context -> current_is_strand_jumped, global_context -> config.space_type);
 
	if(whole_section_matched + sofar_matched > explain_context -> best_matching_bases|| (whole_section_matched + sofar_matched ==  explain_context -> best_matching_bases &&  explain_context -> best_jump_length > explain_context -> tmp_jump_length))
	{
		explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_end =  explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start + remainder_len;
		explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].event_after_section = NULL;

		explain_context -> best_matching_bases = whole_section_matched + sofar_matched ;
		explain_context -> front_search_confirmed_sections = explain_context -> tmp_search_sections +1;
		explain_context -> best_jump_length = explain_context -> tmp_jump_length;
		memcpy(explain_context -> front_search_junctions, explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
	}
}

// read_tail_abs_offset is actually the offset of the base next to the last base in read tail.
// read_tail_pos is the FIRST UNWANTED BASE, after the read.
void search_events_to_back(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context, char * read_text , char * qual_text, unsigned int read_tail_abs_offset, short read_tail_pos, short sofar_matched)
{
	short tested_read_pos;

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;

	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

	int event_search_method;
	if(global_context -> config.do_fusion_detection)
		event_search_method = EVENT_SEARCH_BY_BOTH_SIDES;
	else
		event_search_method = EVENT_SEARCH_BY_LARGE_SIDE;


	// minimum perfect section length is 1
	// tested_read_pos is the first WANTED BASE in section.
	if(MAX_EVENTS_IN_READ - 1> explain_context -> tmp_search_sections)
		for(tested_read_pos = read_tail_pos - 1 - 15 ; tested_read_pos >=0;tested_read_pos --)
		{
			int xk1, matched_bases_to_site;
			int jump_penalty = 0;
			chromosome_event_t *site_events[MAX_EVENT_ENTRIES_PER_SITE];

			int potential_event_pos;

			if(explain_context -> current_is_strand_jumped)
				potential_event_pos = read_tail_abs_offset + ( read_tail_pos - tested_read_pos);
			else
				potential_event_pos = read_tail_abs_offset - ( read_tail_pos - tested_read_pos);
	

			int site_events_no = search_event(global_context, event_table , event_space , potential_event_pos, event_search_method , CHRO_EVENT_TYPE_INDEL | CHRO_EVENT_TYPE_JUNCTION | CHRO_EVENT_TYPE_FUSION , site_events);
			//if(explain_context -> pair_number==2074) printf("BF OFFSET=%d; REDGE=%u; FOUND=%d\n", tested_read_pos, potential_event_pos, site_events_no);


			/*if(memcmp(explain_context->read_name, "HKJMKOB02G7RDV",14) == 0)
			{
				printf("FOUND THE EVENT BACK:%d at %u\n", site_events_no, potential_event_pos);
				if(site_events_no)
					printf("EVENT0_type = %d\n", site_events[0]->event_type);
			}*/

			if(!site_events_no)continue;

			unsigned int tested_chro_begin;
			if(explain_context -> current_is_strand_jumped)
				tested_chro_begin = read_tail_abs_offset + 1;
			else
				tested_chro_begin = read_tail_abs_offset - (read_tail_pos - tested_read_pos);

			matched_bases_to_site = match_chro(read_text + tested_read_pos, value_index, tested_chro_begin , read_tail_pos - tested_read_pos, explain_context -> current_is_strand_jumped, global_context -> config.space_type);

			//if(memcmp(explain_context->read_name, "HKJMKOB02G7RDV",14) == 0)
			//	printf("JUMP?%d > %d\n", (1+matched_bases_to_site)*10000 / (read_tail_pos - tested_read_pos) , 9000);

			if((1+matched_bases_to_site)*10000/(read_tail_pos - tested_read_pos) > 9000)
				for(xk1 = 0; xk1 < site_events_no ; xk1++)
				{
					chromosome_event_t * tested_event = site_events[xk1];
					
					// note that read_tail_pos is the first unwanted base.
					int new_read_tail_pos = tested_read_pos;
					if(tested_event->event_type == CHRO_EVENT_TYPE_INDEL) new_read_tail_pos +=  min(0, tested_event -> indel_length);
					// note that read_tail_abs_offset is the first unwanted base.
					unsigned int new_read_tail_abs_offset;

					if(global_context -> config.do_fusion_detection && tested_event->event_type == CHRO_EVENT_TYPE_FUSION)
					{
						new_read_tail_abs_offset = (potential_event_pos == tested_event -> event_small_side)? tested_event -> event_large_side : tested_event -> event_small_side;
						if(tested_event->is_strand_jumped + explain_context -> current_is_strand_jumped == 1)
							new_read_tail_abs_offset--;
						else
							new_read_tail_abs_offset++;
					}
					else
						new_read_tail_abs_offset = tested_event -> event_small_side + 1;

					//int is_ambiguous = tested_event -> is_ambiguous;

					if(new_read_tail_pos>0)// && (new_read_tail_pos>8 || !is_ambiguous))
					{
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start = tested_read_pos;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].event_after_section = tested_event;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_connected_to_large_side = (potential_event_pos == tested_event -> event_small_side);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].read_pos_end = tested_read_pos + min(0, tested_event->indel_length);
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].abs_offset_for_start = new_read_tail_abs_offset; 
						explain_context -> tmp_jump_length += (tested_event->event_large_side - tested_event->event_small_side);

						//if(tested_event->event_type == CHRO_EVENT_TYPE_FUSION) jump_penalty = 1;
						//else if(tested_event->event_type == CHRO_EVENT_TYPE_JUNCTION) jump_penalty = 1;

						int current_is_jumped = explain_context -> current_is_strand_jumped ;
						if(tested_event -> event_type == CHRO_EVENT_TYPE_FUSION && tested_event -> is_strand_jumped)
							explain_context -> current_is_strand_jumped = !explain_context -> current_is_strand_jumped;
						explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections + 1].is_strand_jumped = explain_context -> current_is_strand_jumped;

						explain_context -> tmp_search_sections ++;

						search_events_to_back(global_context, thread_context, explain_context, read_text , qual_text, new_read_tail_abs_offset , new_read_tail_pos, sofar_matched + matched_bases_to_site - jump_penalty);
						explain_context -> tmp_search_sections --;

						explain_context -> current_is_strand_jumped = current_is_jumped;
						explain_context -> tmp_jump_length -= (tested_event->event_large_side - tested_event->event_small_side);
					
					}
					//if(global_context ->config.limited_tree_scan) break;
				}
			if(global_context ->config.limited_tree_scan && explain_context -> full_read_len <= EXON_LONG_READ_LENGTH) break;
		} 

	int whole_section_matched = match_chro(read_text , value_index, read_tail_abs_offset - (explain_context -> current_is_strand_jumped?-1:read_tail_pos), read_tail_pos , explain_context -> current_is_strand_jumped, global_context -> config.space_type);
 
	if(whole_section_matched + sofar_matched > explain_context -> best_matching_bases || (whole_section_matched + sofar_matched == explain_context -> best_matching_bases && explain_context -> best_jump_length > explain_context -> tmp_jump_length))
	{
		explain_context -> tmp_search_junctions[explain_context -> tmp_search_sections].read_pos_start =  0;
		explain_context -> best_matching_bases = whole_section_matched + sofar_matched ;
		explain_context -> back_search_confirmed_sections = explain_context -> tmp_search_sections +1;
		explain_context -> best_jump_length = explain_context -> tmp_jump_length;
		memcpy(explain_context -> back_search_junctions, explain_context -> tmp_search_junctions , sizeof(perfect_section_in_read_t) * (explain_context -> tmp_search_sections +1)); 
	} 
}

int init_junction_tables(global_context_t * context)
{
	return 0;
}

int destroy_junction_tables(global_context_t * context)
{
	return 0;
}
int init_junction_thread_contexts(global_context_t * global_context, thread_context_t * thread_context, int task)
{
    return 0;
}
int finalise_junction_thread(global_context_t * global_context, thread_context_t * thread_context, int task)
{
    return 0;
}


void insert_big_margin_record(global_context_t * global_context , unsigned char * big_margin_record, unsigned char votes, short read_pos_start, short read_pos_end, int read_len, int is_negative)
{
	unsigned char read_pos_start_2 = (is_negative?read_len -read_pos_end:read_pos_start) ;
	unsigned char read_pos_end_2 = (is_negative?read_len -read_pos_start:read_pos_end);
	assert(votes>0);

	if(read_len>255)
	{
		read_pos_start_2>>=2;
		read_pos_end_2>>=2;
	}

	int xk1;
	for(xk1=0; xk1< global_context->config.big_margin_record_size / 3; xk1++)
	{
		if( votes >= big_margin_record[xk1*3])
			break;
	}
	if(xk1< global_context->config.big_margin_record_size / 3)
	{
		int xk2;
		for(xk2 = global_context->config.big_margin_record_size-4; xk2 >= xk1*3; xk2--)
			big_margin_record[xk2 + 3] = big_margin_record[xk2];
		big_margin_record[xk1*3+0] = votes;
		big_margin_record[xk1*3+1] = read_pos_start_2;
		big_margin_record[xk1*3+2] = read_pos_end_2;
	}
}

//#define voting_anchor_number 3
void set_zero_votes(global_context_t * global_context, int pair_number, int is_second_read , int best_read_id)
{
	if(best_read_id >= global_context->config.multi_best_reads) return;
	_global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id)->selected_votes = 0;
}


void make_128bit_score(unsigned long long int * score_H, unsigned int * score_L, int is_paired_end, short Vote_Anchor_Major, short Vote_Anchor_Minor , short Vote_Second_Major, short Vote_Second_Minor, short Span , short HammingMatch,short Quality, unsigned int TLEM, int Intron)
{
	( * score_H) = 0LLU;
	( * score_L) = 0;

	( * score_H) += (is_paired_end?1LLU:0LLU)<<63; 
	( * score_H) += (1LLU*Vote_Anchor_Major&63)<<57;
	( * score_H) += (1LLU*Vote_Anchor_Minor&31)<<52;
	( * score_H) += (1LLU*Vote_Second_Major&63)<<46;
	( * score_H) += (1LLU*Vote_Second_Minor&31)<<41;

	( * score_H) += (1LLU*Span & 0xfff) << 29;
	( * score_H) += (1LLU*HammingMatch & 0xfff) << 17;
	( * score_H) += (1LLU*Quality & 0x1ff) << 8;
	( * score_H) +=  0xff & (TLEM >> 12);


	( * score_L) += (TLEM & 0xfff) << 20;
	( * score_L) += (Intron & 0xfffff);
}

int process_voting_junction(global_context_t * global_context, thread_context_t * thread_context, int pair_number, gene_vote_t * vote_1, gene_vote_t * vote_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, int read_len_1, int read_len_2, int is_negative_strand)
{
	int i, j, kx1;
	int voting_anchor_number = global_context -> input_reads.is_paired_end_reads?10:global_context -> config.multi_best_reads;

	// each read nominates at most five anchors
	// the base combination of the two anchors is selected.

	select_junction_record_t read_1_anchors[voting_anchor_number];
	select_junction_record_t read_2_anchors[voting_anchor_number];
	int used_anchors_1=0, used_anchors_2=0, is_anchor_1_breakeven = 0, is_anchor_2_breakeven = 0;
	memset(read_1_anchors, 0, sizeof(select_junction_record_t)*voting_anchor_number);
	memset(read_2_anchors, 0, sizeof(select_junction_record_t)*voting_anchor_number);

	int is_second_read;
	int is_junction_found = 0;
	int all_max_votes =  vote_1->max_vote;
	if(global_context -> input_reads.is_paired_end_reads)
		all_max_votes = max(vote_2->max_vote, all_max_votes);


	if(all_max_votes<global_context-> config.minimum_subread_for_first_read)
		return 0;

	for(is_second_read = 0; is_second_read < 1+global_context -> input_reads.is_paired_end_reads; is_second_read++)
	{
		gene_vote_t * current_vote = is_second_read?vote_2:vote_1;
		int current_max_votes = current_vote -> max_vote;
		int total_used_anchors;
		select_junction_record_t * current_anchors = is_second_read?read_2_anchors:read_1_anchors;

		int curr_read_len = is_second_read?read_len_2:read_len_1;
		char * curr_read_text = is_second_read?read_text_2:read_text_1;
		gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

		// put main_piece to anchors.
		for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		{
			for (j=0; j< current_vote->items[i]; j++)
			{
				if(current_vote -> votes[i][j] >=current_max_votes)
				{

					int target_addr = 0;
					int is_break_even ;
					int hamming_match = 0, quality_score = 0;

					if(global_context -> config.use_hamming_distance_break_ties)
						hamming_match = match_chro_indel(curr_read_text, value_index , current_vote -> pos[i][j], curr_read_len, 0, global_context -> config.space_type, global_context -> config.max_indel_length, current_vote -> indel_recorder[i][j], global_context -> config.total_subreads); 
					if(global_context -> config.use_quality_score_break_ties)
						quality_score = max(0,min(512,current_vote -> quality[i][j] / current_vote -> votes[i][j]-200));

					//printf("Q=%d\n", current_vote -> quality[i][j]);

					int main_piece_indels = 0;

					if(curr_read_len > EXON_LONG_READ_LENGTH){
						for(kx1=0; kx1<MAX_INDEL_SECTIONS; kx1++)
						{
							if(!current_vote -> indel_recorder[i][j][kx1*3]) break;
							main_piece_indels += (current_vote -> indel_recorder[i][j][kx1*3+2]);
						}
					}

					unsigned int test_score_L = 0;
					unsigned long long int test_score_H = 0;
					//int test_score = 20000000* current_vote -> votes[i][j] + this_extra_scores + (current_vote -> coverage_end[i][j] - current_vote -> coverage_start[i][j]) - 100 * (main_piece_indels); 
					
					make_128bit_score(&test_score_H, &test_score_L, 0, current_vote -> votes[i][j], 0, 0, 0, (current_vote -> coverage_end[i][j] - current_vote -> coverage_start[i][j]) , hamming_match, quality_score, 0, 0);

					for(target_addr =0; target_addr<voting_anchor_number; target_addr++)
						if((current_anchors[target_addr].Score_H < test_score_H || (current_anchors[target_addr].Score_H == test_score_H && current_anchors[target_addr].Score_L < test_score_L ))|| ( current_vote -> pos[i][j] < current_anchors[target_addr].piece_main_abs_offset && current_anchors[target_addr].Score_H == test_score_H && current_anchors[target_addr].Score_L == test_score_L)) break;

					is_break_even = 0;
					if(current_anchors[0].Score_H == test_score_H && current_anchors[0].Score_L == test_score_L)
						is_break_even = 1;
					else if(current_anchors[0].Score_H < test_score_H || (current_anchors[0].Score_H == test_score_H && current_anchors[0].Score_L < test_score_L))
					{
						if(is_second_read) is_anchor_2_breakeven  = 0;
						else		is_anchor_1_breakeven  = 0;
					}

					if(target_addr<voting_anchor_number-1)
						for(kx1=voting_anchor_number-1; kx1>target_addr; kx1--)
							memcpy(current_anchors+kx1, current_anchors+kx1-1, sizeof(select_junction_record_t));

					if(target_addr<voting_anchor_number)
					{
						memset(&current_anchors[target_addr], 0, sizeof(select_junction_record_t));
						current_anchors[target_addr].piece_main_abs_offset = current_vote -> pos[i][j];
						current_anchors[target_addr].piece_main_coverage_start = current_vote -> coverage_start[i][j];
						current_anchors[target_addr].piece_main_coverage_end   = current_vote -> coverage_end[i][j];
						current_anchors[target_addr].piece_main_votes = current_vote -> votes[i][j];
						current_anchors[target_addr].piece_main_indel_record = current_vote -> indel_recorder[i][j] ;
						current_anchors[target_addr].piece_main_indels = main_piece_indels;
						current_anchors[target_addr].piece_main_masks = current_vote -> masks[i][j];
						current_anchors[target_addr].piece_main_read_quality = quality_score;
						current_anchors[target_addr].piece_main_hamming_match = hamming_match;

						if(global_context -> config.use_hamming_distance_in_exon)
						{
							int found_indels , found_inde_pos;
							
							int matchingness_count = match_indel_chro_to_front(curr_read_text, value_index,  current_vote -> pos[i][j] , curr_read_len, &found_indels, &found_inde_pos, global_context -> config.max_indel_length, 0);

							if(matchingness_count*1000 >= curr_read_len*800)
							{
								current_anchors[target_addr].piece_main_coverage_start = 1;
								current_anchors[target_addr].piece_main_coverage_end = curr_read_len-1;
							}

						}
						current_anchors[target_addr].Score_H = test_score_H;
						current_anchors[target_addr].Score_L = test_score_L;
					}
					if(is_break_even)
					{
						if(is_second_read) is_anchor_2_breakeven  = 1;
						else		is_anchor_1_breakeven  = 1;
					}
				}
				if(current_vote -> votes[i][j] >=current_max_votes-2 && (global_context->config.do_big_margin_filtering_for_junctions || global_context->config.do_big_margin_filtering_for_reads || global_context->config.do_big_margin_reporting))
					insert_big_margin_record(global_context, _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read) ,current_vote -> votes[i][j],  current_vote -> coverage_start[i][j], current_vote -> coverage_end[i][j], is_second_read?read_len_2:read_len_1, is_negative_strand);
			}
		}

		for(kx1=0; kx1<voting_anchor_number; kx1++)
			if(!current_anchors[kx1].piece_main_votes)break;
		total_used_anchors = kx1;

		if(is_second_read)
			used_anchors_2 = total_used_anchors;
		else
			used_anchors_1 = total_used_anchors;

		for(kx1=0; kx1<total_used_anchors; kx1++)
		{
			select_junction_record_t * current_anchor = &current_anchors[kx1];
			//if((current_anchors[kx1].piece_main_coverage_end - current_anchors[kx1].piece_main_coverage_start)*10000 > curr_read_len * 8000)continue;

			if(global_context->config.is_rna_seq_reads || global_context->config.do_fusion_detection)
			{
				unsigned int max_score_L = current_anchor ->Score_L;
				unsigned long long int max_score_H = current_anchor ->Score_H;

				for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
					for (j=0; j< current_vote->items[i]; j++)
					{
						if(current_vote -> pos[i][j] == current_anchor->piece_main_abs_offset) continue;	// myself
						if(current_vote -> votes[i][j] > current_anchor->piece_main_votes) continue;
						if(current_vote -> votes[i][j] == current_anchor->piece_main_votes  && current_vote -> pos[i][j] > current_anchor->piece_main_abs_offset) continue;

						long long int dist = current_vote -> pos[i][j];
						dist -= current_anchor->piece_main_abs_offset;
						int is_strand_jumped = (current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND)!=(current_vote -> masks[i][j] & IS_NEGATIVE_STRAND);

						if(!global_context->config.do_fusion_detection)
						{	// if it is junction detection, then remove long-distance halves and wrongly ordered halves.
							assert(!is_strand_jumped);
							if(abs(dist)> global_context->config.maximum_intron_length) continue; 
							if(current_anchors[kx1].piece_main_coverage_start == current_vote -> coverage_start[i][j])continue;
							if(current_anchors[kx1].piece_main_coverage_end == current_vote -> coverage_end[i][j])continue;

							if(current_anchors[kx1].piece_main_coverage_start > current_vote -> coverage_start[i][j])
							{
								if(current_anchor->piece_main_abs_offset < current_vote -> pos[i][j])continue;
							}
							else
							{
								if(current_anchor->piece_main_abs_offset > current_vote -> pos[i][j])continue;
							}
						}

						int minor_hamming_match = 0;
						if(global_context -> config.use_hamming_distance_break_ties)
							minor_hamming_match = match_chro_indel(curr_read_text, value_index , current_vote -> pos[i][j], curr_read_len, 0, global_context -> config.space_type, global_context -> config.max_indel_length, current_vote -> indel_recorder[i][j], global_context -> config.total_subreads); 

						int minor_read_quality = 0;
						if(global_context -> config.use_quality_score_break_ties)
							minor_read_quality = min(1000, current_vote -> quality[i][j] / current_vote -> votes[i][j]);

						unsigned long long int new_score_H = 0;
						unsigned int new_score_L = 0 ;
						make_128bit_score(&new_score_H, &new_score_L, 0, current_anchor->piece_main_votes , current_vote -> votes[i][j], 0, 0, (current_anchor->piece_main_coverage_end -current_anchor->piece_main_coverage_start) + (current_vote -> coverage_end[i][j] - current_vote -> coverage_start[i][j]) , current_anchors[kx1].piece_main_hamming_match + minor_hamming_match, current_anchors[kx1].piece_main_read_quality + minor_read_quality ,  0 , 1024*1024-1-abs(dist));

						//new_score = current_anchors[kx1].piece_main_extra_scores +  max(500000-abs(dist),0) + current_anchor -> piece_main_votes * 20000000 + current_vote -> votes[i][j] * 20000000 + (current_anchor->piece_main_coverage_end - current_anchor->piece_main_coverage_start) - 100 * (current_anchors[kx1].piece_main_indels);

						if(new_score_H > max_score_H||(new_score_H == max_score_H && new_score_L> max_score_L))
						{
							int final_split_point, is_GT_AG_donors, is_donor_found;
							int donors_found_score;
							int minor_indel_offset=0;

							if(is_strand_jumped)
							{

								// both guess_start and guess_end have to be translated to "reversed" read manner.

								int minor_cover_end_as_reversed = (current_vote -> masks[i][j] & IS_NEGATIVE_STRAND)? current_vote -> coverage_end[i][j]:(curr_read_len - current_vote -> coverage_start[i][j]);
								int minor_cover_start_as_reversed = (current_vote -> masks[i][j] & IS_NEGATIVE_STRAND)? current_vote -> coverage_start[i][j]:(curr_read_len - current_vote -> coverage_end[i][j]);
								int main_cover_end_as_reversed = (current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND)?current_anchors[kx1].piece_main_coverage_end:(curr_read_len - current_anchors[kx1].piece_main_coverage_start);
								int main_cover_start_as_reversed = (current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND)?current_anchors[kx1].piece_main_coverage_start:(curr_read_len - current_anchors[kx1].piece_main_coverage_end);

								// no long overlap
								int overlapped ;
								if(main_cover_start_as_reversed > minor_cover_start_as_reversed)
									overlapped = minor_cover_end_as_reversed - main_cover_start_as_reversed;
								else
									overlapped = main_cover_end_as_reversed - minor_cover_start_as_reversed;

								if(overlapped > 14) continue;

								int guess_start_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed)?
											 (minor_cover_end_as_reversed - 15): (main_cover_end_as_reversed - 15);

								int guess_end_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed)?
											 (main_cover_start_as_reversed + 15): (minor_cover_start_as_reversed + 15);

								int is_left_half_negative = 0 != ((current_anchor->piece_main_abs_offset>current_vote -> pos[i][j]?current_vote -> masks[i][j]:current_anchors[kx1].piece_main_masks)&IS_NEGATIVE_STRAND); 
								int is_right_half_negative = !is_left_half_negative;

								int is_left_on_left_as_reversed = (main_cover_start_as_reversed > minor_cover_start_as_reversed) + (current_anchor->piece_main_abs_offset > current_vote -> pos[i][j]) !=1;

								unsigned int left_half_abs_offset = min(current_vote -> pos[i][j],current_anchor->piece_main_abs_offset);
								unsigned int right_half_abs_offset = max(current_vote -> pos[i][j],current_anchor->piece_main_abs_offset);

								donors_found_score = donor_jumped_score(global_context, thread_context, left_half_abs_offset, right_half_abs_offset , max(0, guess_start_as_reversed) , min( guess_end_as_reversed, curr_read_len),  curr_read_text,  curr_read_len, is_left_half_negative, is_right_half_negative, is_left_on_left_as_reversed , is_second_read,  & final_split_point, & is_GT_AG_donors, & is_donor_found);
							}
							else
							{

								char * chro_name_left, *chro_name_right;
								unsigned int chro_pos_left,chro_pos_right;
								// no long overlap
								int overlapped ;
								if(current_anchors[kx1].piece_main_coverage_start > current_vote -> coverage_start[i][j])
									overlapped = current_vote -> coverage_end[i][j] - current_anchors[kx1].piece_main_coverage_start;
								else
									overlapped = current_anchors[kx1].piece_main_coverage_end - current_vote -> coverage_start[i][j];

								if(overlapped > 14) continue;
								if(abs(dist)<6) continue;
								locate_gene_position( current_anchor->piece_main_abs_offset , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
								locate_gene_position( current_vote -> pos[i][j] , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);
								if(chro_name_right!=chro_name_left) continue;
	
								int guess_start = (current_anchors[kx1].piece_main_coverage_start > current_vote -> coverage_start[i][j])?
											 (current_vote -> coverage_end[i][j] - 15): (current_anchors[kx1].piece_main_coverage_end - 15);

								int guess_end = (current_anchors[kx1].piece_main_coverage_start < current_vote -> coverage_start[i][j])?
											 (current_vote -> coverage_start[i][j] + 15): (current_anchors[kx1].piece_main_coverage_start + 15);

								if(global_context -> config.do_fusion_detection && !(current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND))
									// if for fusion, the current read must have been reversed.
									// hence, it is now changed to "main half" view.
									reverse_read(curr_read_text, curr_read_len, global_context -> config.space_type);

								int normally_arranged = 1!=(current_anchors[kx1].piece_main_coverage_start > current_vote -> coverage_start[i][j]) + (current_anchor->piece_main_abs_offset > current_vote -> pos[i][j]);
								int left_indel_offset=0,  right_indel_offset=0;

								int kx2;
								if(curr_read_len > EXON_LONG_READ_LENGTH){
									for(kx2=0; kx2<MAX_INDEL_SECTIONS; kx2++)
									{
										if(!current_vote -> indel_recorder[i][j][kx2*3]) break;
										minor_indel_offset += (current_vote -> indel_recorder[i][j][kx2*3+2]);
									}
									if(current_anchor->piece_main_abs_offset<  current_vote -> pos[i][j])
									{
										left_indel_offset=current_anchor->piece_main_indels;
										right_indel_offset=minor_indel_offset;
									}
									else
									{
										right_indel_offset=current_anchor->piece_main_indels;
										left_indel_offset=minor_indel_offset;

									}


									// the section having a smaller coordinate will have indel_offset !=0
									// the section having a larger coordiname MUST HAVE indel_offset == 0
									right_indel_offset=0;
								}

								donors_found_score = donor_score(global_context, thread_context, min(current_anchor->piece_main_abs_offset, current_vote -> pos[i][j]),max(current_anchor->piece_main_abs_offset, current_vote -> pos[i][j]), left_indel_offset, right_indel_offset, normally_arranged , max(0, guess_start) , min( guess_end, curr_read_len),  curr_read_text,  curr_read_len, is_second_read,  & final_split_point, & is_GT_AG_donors, & is_donor_found);

								if(global_context -> config.do_fusion_detection && !(current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND))
									// changed back.
									reverse_read(curr_read_text, curr_read_len, global_context -> config.space_type);

							}


//printf("MINORV=%d\tDONOR_FOUND=%d\n", current_vote -> votes[i][j], donors_found_score);


							if(donors_found_score)
							{
								if(global_context -> config.do_fusion_detection && (!current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND) && !is_strand_jumped)
									final_split_point = curr_read_len - final_split_point;

								current_anchors[kx1].piece_minor_abs_offset =  current_vote -> pos[i][j];
								current_anchors[kx1].piece_minor_votes = current_vote -> votes[i][j];
								current_anchors[kx1].piece_minor_coverage_start = current_vote -> coverage_start[i][j];
								current_anchors[kx1].piece_minor_coverage_end = current_vote -> coverage_end[i][j];
								current_anchors[kx1].piece_minor_hamming_match = minor_hamming_match;
								current_anchors[kx1].piece_minor_read_quality = minor_read_quality;
								current_anchors[kx1].piece_minor_indel_offset = minor_indel_offset;
								current_anchors[kx1].intron_length = abs(dist);
								current_anchors[kx1].Score_H = new_score_H;
								current_anchors[kx1].Score_L = new_score_L;
								current_anchors[kx1].split_point = final_split_point;
								current_anchors[kx1].is_GT_AG_donors = is_GT_AG_donors;
								current_anchors[kx1].is_donor_found = is_donor_found;
								if(!is_donor_found)is_junction_found = is_donor_found;
								current_anchors[kx1].is_strand_jumped = is_strand_jumped ;
								max_score_H = new_score_H;
								max_score_L = new_score_L;
							}
						}
					}
			}
			if(current_anchors[kx1].is_strand_jumped)
			{
				// If "is_strand_jumped" is true, all coordinates so far are on the best voted strands (must be differnet strands, namely they're very likely to be overlapped). 
				current_anchors[kx1].piece_minor_coverage_start = curr_read_len - current_anchors[kx1].piece_minor_coverage_end;
				current_anchors[kx1].piece_minor_coverage_end = curr_read_len - current_anchors[kx1].piece_minor_coverage_start;

				// Split_point is now the "negative strand read" view. It has to be changed to "main piece" view
				current_anchors[kx1].split_point = (current_anchors[kx1].piece_main_masks & IS_NEGATIVE_STRAND)?current_anchors[kx1].split_point:(curr_read_len-current_anchors[kx1].split_point);
			}
		}
	}
	
	int is_paired_end_selected = (global_context -> input_reads.is_paired_end_reads && is_result_in_PE( _global_retrieve_alignment_ptr(global_context, pair_number, 0, 0) ));
	int best_read_id_r1 ;
	int best_read_id_r2 =0;
	for(best_read_id_r1=0; best_read_id_r1<global_context->config.multi_best_reads; best_read_id_r1++)
		if(_global_retrieve_alignment_ptr(global_context, pair_number, 0, best_read_id_r1)->selected_votes<1)break;


	if(global_context -> input_reads.is_paired_end_reads)
	{

		for(best_read_id_r2=0; best_read_id_r2<global_context->config.multi_best_reads; best_read_id_r2++)
			if(_global_retrieve_alignment_ptr(global_context, pair_number, 1, best_read_id_r2)->selected_votes<1)break;
	
		for(i=0; i<used_anchors_1; i++)
			for(j=0; j<used_anchors_2; j++)
			{
				long long int dist;
				//int all_votes = read_1_anchors[i].piece_main_votes + read_1_anchors[i].piece_minor_votes + read_2_anchors[j].piece_main_votes + read_2_anchors[j].piece_minor_votes;

				dist = read_1_anchors[i].piece_main_abs_offset;
				dist -= read_2_anchors[j].piece_main_abs_offset;

				if(read_1_anchors[i].piece_main_abs_offset > read_2_anchors[j].piece_main_abs_offset) dist += read_len_1;
				else dist -= read_len_2;

				// the two ends of a segment must conform to the order. 

				unsigned long long int new_score_H = 0;
				unsigned int new_score_L = 0;

				int SUM_COVERAGE = read_1_anchors[i].piece_minor_coverage_end - read_1_anchors[i].piece_minor_coverage_start +
						   read_2_anchors[j].piece_minor_coverage_end - read_2_anchors[j].piece_minor_coverage_start +
						   read_1_anchors[i].piece_main_coverage_end - read_1_anchors[i].piece_main_coverage_start +
						   read_2_anchors[j].piece_main_coverage_end - read_2_anchors[j].piece_main_coverage_start ;

				int SUM_HAMMING =  read_1_anchors[i].piece_main_hamming_match +
						   read_1_anchors[i].piece_main_hamming_match +
						   read_2_anchors[j].piece_minor_hamming_match +
						   read_2_anchors[j].piece_minor_hamming_match ;

				int SUM_QUAL    =  read_1_anchors[i].piece_main_read_quality +
						   read_1_anchors[i].piece_main_read_quality +
						   read_2_anchors[j].piece_minor_read_quality +
						   read_2_anchors[j].piece_minor_read_quality ;

				int SUM_OF_INTRONS = 1024*1024-1 - read_1_anchors[i].intron_length - read_2_anchors[j].intron_length;
				int dist_adjust = max(0, 1024*1024-1-abs(dist));


				int anchor_major_votes = (read_1_anchors[i].piece_main_votes > read_2_anchors[j].piece_main_votes)? read_1_anchors[i].piece_main_votes :read_2_anchors[j].piece_main_votes;
				int anchor_minor_votes = (read_1_anchors[i].piece_main_votes > read_2_anchors[j].piece_main_votes)? read_1_anchors[i].piece_minor_votes :read_2_anchors[j].piece_minor_votes;
				int second_major_votes = (read_1_anchors[i].piece_main_votes > read_2_anchors[j].piece_main_votes)? read_2_anchors[j].piece_main_votes :read_1_anchors[i].piece_main_votes;
				int second_minor_votes = (read_1_anchors[i].piece_main_votes > read_2_anchors[j].piece_main_votes)? read_2_anchors[j].piece_minor_votes :read_1_anchors[i].piece_minor_votes;

				make_128bit_score(&new_score_H, &new_score_L,1, anchor_major_votes, anchor_minor_votes, second_major_votes, second_minor_votes, SUM_COVERAGE , SUM_HAMMING , SUM_QUAL,  dist_adjust , SUM_OF_INTRONS);

				//unsigned int new_score = dist_adjust + read_1_anchors[i].single_score + read_2_anchors[j].single_score - 100*(read_1_anchors[i].piece_main_indels + read_2_anchors[j].piece_main_indels)/2;


				if(global_context->config.is_rna_seq_reads && (read_1_anchors[i].piece_minor_votes || read_2_anchors[j].piece_minor_votes))
				{
					if(((dist < 0 && is_negative_strand) || (dist > 0 && !is_negative_strand)) && !global_context -> config.do_fusion_detection )
						continue;
					if(abs(dist) > global_context->config.maximum_pair_distance + 100000)
						continue;
				}
				else
				{
					if(((dist < 0 && is_negative_strand) || (dist > 0 && !is_negative_strand)) && !global_context -> config.do_fusion_detection )
						continue;
	
					if(abs(dist) > global_context->config.maximum_pair_distance || abs(dist)  < global_context->config.minimum_pair_distance)
						continue;
				}

				alignment_result_t * alignment_1_best = _global_retrieve_alignment_ptr(global_context, pair_number, 0, 0);
				alignment_result_t * alignment_2_best = _global_retrieve_alignment_ptr(global_context, pair_number, 1, 0);

				
				if(new_score_H  > alignment_1_best -> Score_H || (new_score_H == alignment_1_best-> Score_H && new_score_L >= alignment_1_best-> Score_L))
				{
					if(new_score_H > alignment_1_best-> Score_H || new_score_L > alignment_1_best-> Score_L)
					{
						best_read_id_r1 = 0;
						best_read_id_r2 = 0;


						alignment_1_best -> result_flags &= ~CORE_IS_BREAKEVEN; 
						alignment_2_best -> result_flags &= ~CORE_IS_BREAKEVEN; 
					}
					else
					{
						//printf("SET_BE: %d ; S=%16llx+%16llX\n", pair_number, new_score_H, new_score_L);
						alignment_1_best -> result_flags |= CORE_IS_BREAKEVEN; 
						alignment_2_best -> result_flags |= CORE_IS_BREAKEVEN; 
					}

					int r1_used_subreads = max(vote_1 -> all_used_subreads, alignment_1_best->used_subreads_in_vote );
					int r2_used_subreads = max(vote_2 -> all_used_subreads, alignment_2_best->used_subreads_in_vote );

					set_alignment_result(global_context, pair_number, 0, best_read_id_r1, read_1_anchors[i].piece_main_abs_offset, read_1_anchors[i].piece_main_votes , read_1_anchors[i].piece_main_indel_record, read_1_anchors[i].piece_main_coverage_start, read_1_anchors[i].piece_main_coverage_end, 0!=(read_1_anchors[i].piece_main_masks & IS_NEGATIVE_STRAND), read_1_anchors[i].piece_minor_abs_offset, read_1_anchors[i].piece_minor_votes, read_1_anchors[i].piece_minor_coverage_start, read_1_anchors[i].piece_minor_coverage_end, read_1_anchors[i].split_point, read_1_anchors[i].is_strand_jumped, read_1_anchors[i].is_donor_found?read_1_anchors[i].is_GT_AG_donors:-1, r1_used_subreads, vote_1-> noninformative_subreads, read_1_anchors[i].piece_main_indels, read_1_anchors[i].piece_minor_indel_offset);
					set_alignment_result(global_context, pair_number, 1, best_read_id_r2, read_2_anchors[j].piece_main_abs_offset, read_2_anchors[j].piece_main_votes , read_2_anchors[j].piece_main_indel_record, read_2_anchors[j].piece_main_coverage_start, read_2_anchors[j].piece_main_coverage_end, 0!=(read_2_anchors[j].piece_main_masks & IS_NEGATIVE_STRAND), read_2_anchors[j].piece_minor_abs_offset, read_2_anchors[j].piece_minor_votes, read_2_anchors[j].piece_minor_coverage_start, read_2_anchors[j].piece_minor_coverage_end, read_2_anchors[j].split_point, read_2_anchors[j].is_strand_jumped, read_2_anchors[j].is_donor_found?read_2_anchors[j].is_GT_AG_donors:-1, r2_used_subreads, vote_2-> noninformative_subreads, read_2_anchors[j].piece_main_indels, read_2_anchors[j].piece_minor_indel_offset);

					alignment_1_best -> Score_H = new_score_H;
					alignment_1_best -> Score_L = new_score_L;
					alignment_2_best -> Score_H = new_score_H;
					alignment_2_best -> Score_L = new_score_L;
					
					is_paired_end_selected = 1;

					best_read_id_r1 += 1;
					best_read_id_r2 += 1;

					set_zero_votes(global_context, pair_number,0 , best_read_id_r1);
					set_zero_votes(global_context, pair_number,1 , best_read_id_r2);
				}
			}
	}

	if(!is_paired_end_selected)
	{
		alignment_result_t * alignment_1_best = _global_retrieve_alignment_ptr(global_context, pair_number, 0, 0);
		for(i=0; i<used_anchors_1; i++)
		{
			if((read_1_anchors[i].Score_H > alignment_1_best -> Score_H) || (read_1_anchors[i].Score_H == alignment_1_best -> Score_H && read_1_anchors[i].Score_L >= alignment_1_best -> Score_L))
			{
				if(read_1_anchors[i].Score_H > alignment_1_best -> Score_H ||  read_1_anchors[i].Score_L > alignment_1_best -> Score_L )
				{
					best_read_id_r1 = 0;

					if(is_anchor_1_breakeven)
						alignment_1_best -> result_flags |= CORE_IS_BREAKEVEN;
					else
						alignment_1_best -> result_flags &= ~CORE_IS_BREAKEVEN;
				}
				else
				{
					if(read_1_anchors[i].piece_main_abs_offset > _global_retrieve_alignment_ptr(global_context, pair_number, 0, 0)->selected_position && global_context->config.multi_best_reads == 1)
						best_read_id_r1 = 0;

					alignment_1_best -> result_flags |= CORE_IS_BREAKEVEN;
				}

				alignment_1_best -> Score_H = read_1_anchors[i].Score_H;
				alignment_1_best -> Score_L = read_1_anchors[i].Score_L;

				//printf("BEST_ID_R1=%d\n",best_read_id_r1);

				// TODO: add result at best_read_id_r1
				alignment_result_t * r1_result = _global_retrieve_alignment_ptr(global_context, pair_number, 0, best_read_id_r1);
				int r1_used_subreads = max(vote_1 -> all_used_subreads, r1_result->used_subreads_in_vote );
				set_alignment_result(global_context, pair_number, 0, best_read_id_r1, read_1_anchors[i].piece_main_abs_offset, read_1_anchors[i].piece_main_votes , read_1_anchors[i].piece_main_indel_record, read_1_anchors[i].piece_main_coverage_start, read_1_anchors[i].piece_main_coverage_end, 0!=(read_1_anchors[i].piece_main_masks & IS_NEGATIVE_STRAND), read_1_anchors[i].piece_minor_abs_offset, read_1_anchors[i].piece_minor_votes, read_1_anchors[i].piece_minor_coverage_start, read_1_anchors[i].piece_minor_coverage_end, read_1_anchors[i].split_point, read_1_anchors[i].is_strand_jumped, read_1_anchors[i].is_donor_found?read_1_anchors[i].is_GT_AG_donors:-1, r1_used_subreads, vote_1-> noninformative_subreads, read_1_anchors[i].piece_main_indels, read_1_anchors[i].piece_minor_indel_offset);

				best_read_id_r1 += 1;
				set_zero_votes(global_context, pair_number,0 , best_read_id_r1);
				/*if(memcmp(read_name_1, "HKJMKOB02G7RDV",14) == 0){

					printf("MAX_VOTES=%d\t\tSTART_POS=%u\t\tMINOR_VOTES=%d\t\tMINOR_POS=%u\n",  read_1_anchors[i].piece_main_votes, read_1_anchors[i].piece_main_abs_offset, read_1_anchors[i].piece_minor_votes, read_1_anchors[i].piece_minor_abs_offset);
					print_votes(vote_1, global_context -> config.index_prefix);
				}*/
			}
		}


		if(global_context -> input_reads.is_paired_end_reads)
		{

			alignment_result_t * alignment_2_best = _global_retrieve_alignment_ptr(global_context, pair_number, 1, 0);
			for(j=0; j<used_anchors_2; j++)
			{
				if(read_2_anchors[j].Score_H > alignment_2_best -> Score_H || (read_2_anchors[j].Score_H == alignment_2_best -> Score_H && read_2_anchors[j].Score_L >= alignment_2_best -> Score_L))
				{
					if(read_2_anchors[j].Score_H > alignment_2_best -> Score_H || read_2_anchors[j].Score_L >  alignment_2_best-> Score_L )
					{
						best_read_id_r2 = 0;

						if(is_anchor_2_breakeven)
							alignment_2_best -> result_flags |= CORE_IS_BREAKEVEN;
						else
							alignment_2_best -> result_flags &= ~CORE_IS_BREAKEVEN;
					}
					else
					{
						if(read_2_anchors[j].piece_main_abs_offset > _global_retrieve_alignment_ptr(global_context, pair_number, 1, 0)->selected_position && global_context->config.multi_best_reads == 1)
							best_read_id_r2 = 0;
						alignment_2_best -> result_flags &= ~CORE_IS_BREAKEVEN;
						//printf("SET_BE_2: %d\n", pair_number);
					}

					alignment_2_best -> Score_H = read_2_anchors[j].Score_H;
					alignment_2_best -> Score_L = read_2_anchors[j].Score_L;
					is_paired_end_selected = 0;

					// TODO: add result at best_read_id_r2
					alignment_result_t * r2_result = _global_retrieve_alignment_ptr(global_context, pair_number, 1, 0);
					int r2_used_subreads = max(vote_2 -> all_used_subreads, r2_result->used_subreads_in_vote );
					set_alignment_result(global_context, pair_number, 1, best_read_id_r2, read_2_anchors[j].piece_main_abs_offset, read_2_anchors[j].piece_main_votes , read_2_anchors[j].piece_main_indel_record, read_2_anchors[j].piece_main_coverage_start, read_2_anchors[j].piece_main_coverage_end, 0!=(read_2_anchors[j].piece_main_masks & IS_NEGATIVE_STRAND), read_2_anchors[j].piece_minor_abs_offset, read_2_anchors[j].piece_minor_votes, read_2_anchors[j].piece_minor_coverage_start, read_2_anchors[j].piece_minor_coverage_end, read_2_anchors[j].split_point, read_2_anchors[j].is_strand_jumped, read_2_anchors[j].is_donor_found?read_2_anchors[j].is_GT_AG_donors:-1, r2_used_subreads, vote_2-> noninformative_subreads,  read_2_anchors[j].piece_main_indels, read_2_anchors[j].piece_minor_indel_offset);

					best_read_id_r2 += 1;
					set_zero_votes(global_context, pair_number,1 , best_read_id_r2);
				}
			}
		}

	}


	alignment_result_t * tmp_result = _global_retrieve_alignment_ptr(global_context, pair_number, 0, 0);
	if(tmp_result->selected_votes <1)
	{
		tmp_result -> used_subreads_in_vote = max(vote_1 -> all_used_subreads, tmp_result -> used_subreads_in_vote );
		tmp_result -> noninformative_subreads_in_vote = max(vote_1 -> noninformative_subreads, tmp_result -> noninformative_subreads_in_vote);
	}


	tmp_result = _global_retrieve_alignment_ptr(global_context, pair_number, 1, 0);
	if(tmp_result->selected_votes <1 && global_context -> input_reads.is_paired_end_reads)
	{

		tmp_result->used_subreads_in_vote = max(vote_2 -> all_used_subreads, tmp_result->used_subreads_in_vote );
		tmp_result->noninformative_subreads_in_vote = max(vote_2 -> noninformative_subreads, tmp_result->noninformative_subreads_in_vote);
	}
	return 0;
}

int explain_read(global_context_t * global_context, thread_context_t * thread_context, int pair_number, int read_len, char * read_name , char *read_text, char *qual_text, int is_second_read, int best_read_id, int is_negative_strand)
{
	explain_context_t explain_context;

	alignment_result_t *current_result = _global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id); 

	if(global_context -> config.do_big_margin_reporting || global_context -> config.do_big_margin_filtering_for_reads)
	{
		int current_repeated_times = is_ambiguous_voting(global_context, pair_number, is_second_read, current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, read_len, (current_result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);
		if(current_repeated_times>1) return 0;
	}
	


	memset(&explain_context,0, sizeof(explain_context_t));

	explain_context.full_read_len = read_len;
	explain_context.full_read_text = read_text;
	explain_context.full_qual_text = qual_text;
	explain_context.read_name = read_name;
	explain_context.is_confirmed_section_negative_strand = is_negative_strand ;
	explain_context.pair_number = pair_number;
	explain_context.is_second_read = is_second_read ;
	explain_context.best_read_id = best_read_id;


	unsigned int back_search_tail_position, front_search_start_position;
	unsigned short back_search_read_tail, front_search_read_start;


	back_search_read_tail = min(explain_context.full_read_len , current_result -> confident_coverage_end );//- 5;
	back_search_tail_position = current_result -> selected_position + back_search_read_tail +  current_result -> indels_in_confident_coverage;

	explain_context.tmp_search_junctions[0].read_pos_end = back_search_read_tail;
	explain_context.tmp_search_junctions[0].abs_offset_for_start = back_search_tail_position;

	explain_context.tmp_jump_length = 0;
	explain_context.best_jump_length = 0xffff0000;

	search_events_to_back(global_context, thread_context, &explain_context, read_text , qual_text, back_search_tail_position , back_search_read_tail, 0);

	if(explain_context.back_search_confirmed_sections>0)
	{
		
		short last_section_length = explain_context.back_search_junctions[0].read_pos_end - explain_context.back_search_junctions[0].read_pos_start;
		
		front_search_read_start = explain_context.back_search_junctions[0].read_pos_start; 
		front_search_start_position = explain_context.back_search_junctions[0].abs_offset_for_start - last_section_length;

		int last_sec = explain_context.back_search_confirmed_sections-1;

		current_result -> selected_position = explain_context.back_search_junctions[last_sec].abs_offset_for_start - explain_context.back_search_junctions[last_sec].read_pos_end + explain_context.back_search_junctions[last_sec].read_pos_start;
 
	}
	else
	{
		front_search_read_start = current_result -> confident_coverage_start + 5;
		front_search_start_position = current_result -> selected_position + front_search_read_start;
	}


	// clean the temporary results
	explain_context.tmp_search_sections = 0;
	explain_context.best_matching_bases = 0;
	memset(explain_context.tmp_search_junctions, 0, sizeof(perfect_section_in_read_t ) * MAX_EVENTS_IN_READ);

	explain_context.tmp_search_junctions[0].read_pos_start = front_search_read_start;
	explain_context.tmp_search_junctions[0].abs_offset_for_start = front_search_start_position;
	explain_context.tmp_jump_length = 0;
	explain_context.best_jump_length = 0xffff0000;
	search_events_to_front(global_context, thread_context, &explain_context, read_text + front_search_read_start, qual_text + front_search_read_start, front_search_start_position,read_len - front_search_read_start , 0);

	// calc
	finalise_explain_CIGAR(global_context, thread_context, &explain_context);

	return 0;
}

int find_soft_clipping(global_context_t * global_context,  thread_context_t * thread_context, gene_value_index_t * current_value_index, char * read_text, unsigned int mapped_pos, int test_len,  int search_to_tail)
{
	#define SOFT_CLIPPING_WINDOW_SIZE 6
	#define SOFT_CLIPPING_MAX_ERROR   1

	char window_matched[SOFT_CLIPPING_WINDOW_SIZE];
	int x0,x1,x2;

	memset(window_matched, 0 , SOFT_CLIPPING_WINDOW_SIZE);

	for(x0=0;x0 < test_len; x0++)
	{

		if(search_to_tail) x1 = test_len -1 -x0;
		else	x1=x0;
		char ref_value = gvindex_get(current_value_index, mapped_pos + x1);
		int sum_matched=0;
		for(x2 = SOFT_CLIPPING_WINDOW_SIZE - 1; x2 > 0; x2--)
		{
			window_matched[x2] = window_matched[x2-1];
			sum_matched += window_matched[x2];
		}
		window_matched[0] = (ref_value == read_text[x1]);
		sum_matched += window_matched[0];

		// find the first matched base, such that the matched bases >= SOFT_CLIPPING_WINDOW_SIZE - SOFT_CLIPPING_MAX_ERROR if this base is added into the window.
		if(window_matched[0])
		{
			if(sum_matched > SOFT_CLIPPING_WINDOW_SIZE - SOFT_CLIPPING_MAX_ERROR)
			{
				return max(0 , x0 - SOFT_CLIPPING_WINDOW_SIZE);
			}
		}
		
	}
	return 0;
}

// read_head_abs_offset is the first WANTED base in read.
// If the first section in read is reversed, read_head_abs_offset is the LAST WANTED bases in this section. (the abs offset of the first base in the section is actually larger than read_head_abs_offset)
int final_CIGAR_quality(global_context_t * global_context, thread_context_t * thread_context, char * read_text, char * qual_text, int read_len, char * cigar_string, unsigned long read_head_abs_offset, int is_read_head_reversed, int * mismatched_bases)
{
	int cigar_cursor = 0;
	int read_cursor = 0;
	unsigned int current_perfect_section_abs = read_head_abs_offset;
	int rebuilt_read_len = 0;
	float all_matched_bases = 0;
	gene_value_index_t * current_value_index = thread_context?thread_context->current_value_index:global_context->current_value_index; 
	int current_reversed = is_read_head_reversed;
	int all_perfect_length = 0;
	int is_First_M = 1;
	int head_soft_clipped = -1, tail_soft_clipped = -1;

	unsigned int tmp_int = 0;
	while(1)
	{
		char nch = cigar_string[cigar_cursor++];
		if(!nch)break;
		if(isdigit(nch))
			tmp_int = tmp_int*10+(nch-'0');
		else{
			if(nch == 'M' || nch == 'S')
			{
				char *qual_text_cur;
				if(qual_text[0])qual_text_cur = qual_text+read_cursor;
				else	qual_text_cur = NULL;

				float section_qual = match_base_quality(current_value_index, read_text+read_cursor, current_perfect_section_abs, qual_text_cur, tmp_int, current_reversed, global_context->config.phred_score_format , mismatched_bases, global_context -> config.high_quality_base_threshold);
				all_matched_bases += section_qual;
				rebuilt_read_len += tmp_int;
				all_perfect_length += tmp_int;

				int is_Last_M = (cigar_string[cigar_cursor]==0);

				// find "J" sections if it is the first M
				if(is_First_M && global_context -> config.show_soft_cliping)
				{
					head_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, read_text, current_perfect_section_abs, tmp_int, 0);
					if(head_soft_clipped == tmp_int) head_soft_clipped = 0;
				}
				if(is_Last_M && global_context -> config.show_soft_cliping)
				{
					tail_soft_clipped = find_soft_clipping(global_context, thread_context, current_value_index, read_text + read_cursor, current_perfect_section_abs, tmp_int, 1);
					if(tail_soft_clipped == tmp_int) tail_soft_clipped = 0;
				}
				if(is_Last_M && is_First_M && tail_soft_clipped+head_soft_clipped >= tmp_int-1)
				{
					head_soft_clipped=0;
					tail_soft_clipped=0;
				}
				is_First_M=0;


				read_cursor += tmp_int;

				//move to the NEXT UNWANTED ABS OFFSET. 
				if(current_reversed)
					current_perfect_section_abs --;
				else
					current_perfect_section_abs += tmp_int;


			}
			else if(nch == 'I')
			{
				rebuilt_read_len += tmp_int;
				read_cursor += tmp_int;
				all_matched_bases += tmp_int;
			}
			else if(nch == 'D')
			{
				if(!current_reversed)
					current_perfect_section_abs += tmp_int;
			}
			else if(tolower(nch) == 'n')
			{
				current_perfect_section_abs += tmp_int;
				if(nch == 'n') current_reversed = !current_reversed;
			}
			else if(tolower(nch) == 'b')
			{
				current_perfect_section_abs -= tmp_int;
				if(nch == 'b') current_reversed = !current_reversed;
			}
			tmp_int = 0;
		}
	}

	assert(rebuilt_read_len == read_len);


	if(global_context -> config.show_soft_cliping && (head_soft_clipped>0 || tail_soft_clipped>0))
	{
		char new_cigar_tmp[100];
		is_First_M=1;
		new_cigar_tmp[0]=0;
		cigar_cursor = 0;
		while(1)
		{
			char nch = cigar_string[cigar_cursor++];

			if(!nch)break;
			if(isdigit(nch))
				tmp_int = tmp_int*10+(nch-'0');
			else{
				char cigar_piece [30];
				cigar_piece[0]=0;

				if(nch == 'M')
				{
					char cigar_tiny [11];
					int is_Last_M = (cigar_string[cigar_cursor]==0);
					if(is_First_M && head_soft_clipped>0)
					{
						tmp_int -= head_soft_clipped;
						sprintf(cigar_tiny,"%dS",head_soft_clipped);
						strcat(cigar_piece, cigar_tiny);
					}
					if(is_Last_M && tail_soft_clipped>0)
					{
						tmp_int -= tail_soft_clipped;
					}
					sprintf(cigar_tiny,"%dM",tmp_int);
					strcat(cigar_piece, cigar_tiny);
					if(is_Last_M && tail_soft_clipped>0)
					{
						sprintf(cigar_tiny,"%dS",tail_soft_clipped);
						strcat(cigar_piece, cigar_tiny);
					}
					is_First_M = 0;
				}
				else
				{
					sprintf(cigar_piece, "%u%c", tmp_int, nch);
				}

				strcat(new_cigar_tmp, cigar_piece);
				tmp_int = 0;
			}
		}

		strcpy(cigar_string, new_cigar_tmp);
	}

	return 100+(int)(all_matched_bases*100/read_len);
}

// this function also adds final_counting_reads in chromosome_events.
int finalise_explain_CIGAR(global_context_t * global_context, thread_context_t * thread_context, explain_context_t * explain_context)
{
	int xk1;
	char tmp_cigar[100];
	chromosome_event_t * to_be_supported [20];
	short flanking_size_left[20], flanking_size_right[20];
	int to_be_supported_count = 0;
	int is_junction_read = 0;
	int total_perfect_matched_sections = 0;
	alignment_result_t * result = _global_retrieve_alignment_ptr(global_context, explain_context->pair_number, explain_context->is_second_read, explain_context-> best_read_id); 


	tmp_cigar[0]=0;
	// reverse the back_search results
	for(xk1=0; xk1<explain_context -> back_search_confirmed_sections/2; xk1++)
	{
		perfect_section_in_read_t tmp_exp;
		memcpy(&tmp_exp, &explain_context -> back_search_junctions[xk1], sizeof(perfect_section_in_read_t));
		memcpy(&explain_context -> back_search_junctions[xk1],  &explain_context -> back_search_junctions[explain_context -> back_search_confirmed_sections - xk1 - 1] , sizeof(perfect_section_in_read_t));
		memcpy(&explain_context -> back_search_junctions[explain_context -> back_search_confirmed_sections - xk1 - 1] , &tmp_exp , sizeof(perfect_section_in_read_t));
	} 
	
	// adding indel lengths in read lengths and relocate sections
	// note that the last section in back results has the same strand of the main piece.
	int is_first_section_negative = (result ->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0; 
	for(xk1=0; xk1<explain_context -> back_search_confirmed_sections; xk1++)
	{
		int section_length = explain_context -> back_search_junctions[xk1].read_pos_end - explain_context -> back_search_junctions[xk1].read_pos_start; 
		unsigned int new_start_pos;

		if(explain_context -> back_search_junctions[xk1].is_strand_jumped)
			// the "strand_jumped" section do not need to move
			// however, the "abs_offset_for_start" is actually for the last base in this section.
			// this does not metter if we compare the reversed read to the chromosome.
			// "abs_offset_for_start" is the first UNWANTED base (smaller than the first WANTED base)
			new_start_pos = explain_context -> back_search_junctions[xk1].abs_offset_for_start +1;
		else
			// "abs_offset_for_start" is the first UNWANTED base. By subtracting the length, it becomes the first WANTED base.
			new_start_pos = explain_context -> back_search_junctions[xk1].abs_offset_for_start - section_length;

		explain_context -> back_search_junctions[xk1].abs_offset_for_start = new_start_pos;
		if(explain_context -> back_search_junctions[xk1].event_after_section
			&& explain_context -> back_search_junctions[xk1].event_after_section->is_strand_jumped) is_first_section_negative=!is_first_section_negative;
	}

	// build CIGAR
	int is_cigar_overflow = 0;
	for(xk1 = 0; xk1 < explain_context -> back_search_confirmed_sections + explain_context -> front_search_confirmed_sections -1; xk1++)
	{
		char piece_cigar[20];
		int read_pos_start, read_pos_end;
		perfect_section_in_read_t * current_section, *next_section = NULL;

		int is_front_search = 0;
		if(xk1 >= explain_context -> back_search_confirmed_sections || xk1 == explain_context -> back_search_confirmed_sections -1)
		{
			current_section = &explain_context -> front_search_junctions[xk1 - explain_context -> back_search_confirmed_sections +1];
			if(xk1 - explain_context -> back_search_confirmed_sections +2 < explain_context -> front_search_confirmed_sections)
				next_section = &explain_context -> front_search_junctions[xk1 - explain_context -> back_search_confirmed_sections +2];

			is_front_search = 1;
		}
		else
		{
			current_section = &explain_context -> back_search_junctions[xk1];
			if(xk1+1 <  explain_context -> back_search_confirmed_sections)
				next_section = &explain_context -> back_search_junctions[xk1+1];
		}


		read_pos_start = current_section -> read_pos_start;
		read_pos_end = current_section -> read_pos_end;
		chromosome_event_t *event_after = current_section -> event_after_section;

		sprintf(piece_cigar, "%dM", (read_pos_end - read_pos_start));
		total_perfect_matched_sections += (read_pos_end - read_pos_start);
		flanking_size_left[xk1] = (read_pos_end - read_pos_start);

		if(xk1<explain_context -> back_search_confirmed_sections + explain_context -> front_search_confirmed_sections -2)
			assert(event_after);

		if(xk1>0)
			flanking_size_right[xk1-1] = (read_pos_end - read_pos_start);

		if(event_after)
		{
			if(event_after -> event_type == CHRO_EVENT_TYPE_INDEL)
				sprintf(piece_cigar+strlen(piece_cigar), "%d%c", abs(event_after->indel_length), event_after->indel_length>0?'D':'I');
			else if(event_after -> event_type == CHRO_EVENT_TYPE_JUNCTION||event_after -> event_type == CHRO_EVENT_TYPE_FUSION)
			{
				char jump_mode = current_section -> is_connected_to_large_side?'B':'N';
				if(event_after -> is_strand_jumped) jump_mode = tolower(jump_mode);

				// the distance in CIGAR is the NEXT UNWANTED BASE of piece#1 to the FIRST WANTED BASE in piece#2
				int delta_one ;
				if(current_section -> is_strand_jumped + current_section -> is_connected_to_large_side == 1) delta_one = 1;
				else delta_one = -1;

				// if it is from front_search, the event side points to the first WANTED base of the next section; it should be moved to the last WANTED base the next section if the next section is jumped.
				if(next_section && (event_after -> is_strand_jumped + current_section -> is_strand_jumped==1))
				{
					if(is_front_search)
					{
						if(current_section -> is_connected_to_large_side)
							delta_one += (next_section->read_pos_end - next_section-> read_pos_start - 1);
						else
							delta_one -= (next_section->read_pos_end - next_section-> read_pos_start - 1);
					}
					else
					{
						if(current_section -> is_connected_to_large_side)
							delta_one += (next_section->read_pos_end - next_section-> read_pos_start - 1);
						else
							delta_one -= (next_section->read_pos_end - next_section-> read_pos_start - 1);
					}
				}
				
				sprintf(piece_cigar+strlen(piece_cigar), "%d%c", event_after -> event_large_side - event_after -> event_small_side + delta_one, jump_mode);
				is_junction_read ++;
			}
			to_be_supported[to_be_supported_count++] = event_after;
		}
		strcat(tmp_cigar, piece_cigar);
		if(strlen(tmp_cigar)>80){
			is_cigar_overflow=1;
			break;
		}
	}

	int mismatch_bases = 0, isCigarOK = 0;

	if(is_cigar_overflow) sprintf(tmp_cigar, "%dM",  explain_context -> full_read_len);

	unsigned int final_position = explain_context -> back_search_junctions[0].abs_offset_for_start;
	int final_qual = final_CIGAR_quality(global_context, thread_context, explain_context -> full_read_text, explain_context -> full_qual_text, explain_context -> full_read_len , tmp_cigar, final_position, is_first_section_negative != ((result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0), &mismatch_bases);

	//if(memcmp(explain_context->read_name, "HKJMKOB02G7RDV",14) == 0)printf("POS=%u\tCIGAR=%s\tMM=%d\tQUAL=%d\n", final_position , tmp_cigar, mismatch_bases, final_qual);

	int applied_mismatch = is_junction_read? global_context->config.max_mismatch_junction_reads:global_context->config.max_mismatch_exonic_reads ;
	if(explain_context->full_read_len > EXON_LONG_READ_LENGTH)
		applied_mismatch = ((((explain_context->full_read_len+1)<<16) / 100) * applied_mismatch)>>16;
	if(mismatch_bases <= applied_mismatch)
	{
		int compressed_len;
		result -> final_quality = final_qual;
		result -> selected_position = final_position;
		if(((result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) != is_first_section_negative)
		{
			assert(0);
			result -> cigar_string[0]=0xff;
			compressed_len = cigar2bincigar(tmp_cigar, result -> cigar_string + 1, CORE_MAX_CIGAR_LEN - 1);
		}
		else
			compressed_len = cigar2bincigar(tmp_cigar, result -> cigar_string, CORE_MAX_CIGAR_LEN);

		// commit the change to the chromosome_events
		if(compressed_len>0)
		{
			for(xk1= 0; xk1 < to_be_supported_count; xk1++)
			{
				if(thread_context)
					((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> final_counted_reads_array [ to_be_supported [xk1] -> global_event_id] ++;
				else
					to_be_supported [xk1] -> final_counted_reads ++;
				if(to_be_supported [xk1] -> event_type !=CHRO_EVENT_TYPE_INDEL) 
				{
					short current_event_flanking_left = flanking_size_left[xk1];
					short current_event_flanking_right = flanking_size_right[xk1];
					to_be_supported [xk1] -> junction_flanking_left = max(to_be_supported [xk1] -> junction_flanking_left, current_event_flanking_left); 
					to_be_supported [xk1] -> junction_flanking_right = max(to_be_supported [xk1] -> junction_flanking_right, current_event_flanking_right); 
				}
			}

			result -> result_flags |= CORE_IS_FULLY_EXPLAINED;
			isCigarOK=1;
		}
	}
	//printf("BRRRRRID=%d; MM=%d;   CIGAR=%s<%d>   QUAL=%d\n",  explain_context-> best_read_id, mismatch_bases, tmp_cigar, isCigarOK, final_qual);

	if(!isCigarOK)
	{
		result -> final_quality = final_qual;
		result -> result_flags &= ~CORE_IS_FULLY_EXPLAINED;
		result -> Score_H &= 0x7fffffffffffffffllu;
	}
	return 0;
}




#define ceq(c,t) ((c)[0]==(t)[0] && (c)[1]==(t)[1])
#define c2eq(ch1, ch2, tg1, tg2) ((ceq(ch1, tg1) && ceq(ch2, tg2)) || (ceq(ch1, tg2) && ceq(ch2, tg1)) )

int paired_chars_full_core(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 2;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 2;
	}
	else if ( c2eq(ch1, ch2,"GC","AG") || c2eq(ch1, ch2,"GC","CT") || c2eq(ch1, ch2,"AT","AC") || c2eq(ch1, ch2,"GT","AT"))
	{
		if (is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT")  || ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "GC") || ceq(ch1, "AT") ||ceq(ch1, "GT") || ceq(ch1, "CT")) return 1;
	}
	return 0;
}

int paired_chars_part_core(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse)
		{
			if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		}else
			if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
	}
	return 0;
}

#define  paired_chars paired_chars_part_core


#define is_donor_chars_full(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T') || \
			    ((cc)[0]=='G' && (cc)[1]=='C') || \
			    ((cc)[0]=='A' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') ) 


#define is_donor_chars_part(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T')) 

#define is_donor_chars is_donor_chars_part




int is_ambiguous_voting(global_context_t * global_context, int pair_number, int is_second_read, int max_vote, int max_start,int max_end, int read_len, int is_negative)
{
	int xk1;
	int encounter = 0;

	if(is_negative)
	{
		int tmp = max_start;
		max_start = read_len - max_end;
		max_end = read_len - tmp;
	}

	if(read_len > 255)
	{
		max_start = max_start>>2;
		max_end = max_end>>2;
	}

	unsigned char * big_margin_record = _global_retrieve_big_margin_ptr(global_context,pair_number, is_second_read);

	for(xk1 = 0; xk1 < global_context->config.big_margin_record_size/3 ; xk1++)
	{
		if(!big_margin_record[xk1*3])break;

		if((big_margin_record[xk1*3]) >= max_vote -1)	// actually, max-1
			if(big_margin_record[xk1*3+1] >= max_start - 2 && big_margin_record[xk1*3+2] <= max_end + 1)
				encounter++;

	}
	if(encounter>1) return encounter;
	return 0;
}

#define JUNCTION_CONFIRM_WINDOW 17
// This function implements the same function of donor_score, except that the two halves are from different strands.
// Both halves are forced to positive strand and the split point is found.
// Note that the donor/receptor sides are still expected for distinguishing between Fusion Breaks and Fusion Junctions.

// Note that the read_text is on reversed mode. The guess points are on reversed mode too.
// "Left" and "Right" means the left/right half in the "reversed" read.
int donor_jumped_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int left_virtualHead_abs_offset, unsigned int right_virtualHead_abs_offset, int guess_start, int guess_end,  char * read_text, int read_len, int is_left_half_negative, int is_right_half_negative, int normally_arranged, int is_second_read, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found)
{
	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	// guess_end is the index of the first UNWANTED BASE.
	int most_likely_point_as_reversed = (guess_start+guess_end)/2;

	int selected_real_split_point = -1, selected_junction_strand = -1;
	//char donor_left[2], donor_right[2];
 
	int best_score = -111111;

	int real_split_point_i;
	int real_split_point_numbers = guess_end - guess_start;

	char positive_read[MAX_READ_LENGTH+1];
	strcpy(positive_read, read_text) ;
	reverse_read(positive_read, read_len, global_context->config.space_type);
	
	for(real_split_point_i = 0 ; real_split_point_i < real_split_point_numbers; real_split_point_i++)
	{
		int left_should_match, right_should_match;
		int left_should_not_match, right_should_not_match;
		int real_split_point_as_reversed = (real_split_point_i % 2)?-((real_split_point_i+1)/2):((1+real_split_point_i)/2);
		real_split_point_as_reversed += most_likely_point_as_reversed;

		if(real_split_point_as_reversed > read_len-JUNCTION_CONFIRM_WINDOW)continue;
		if(real_split_point_as_reversed < JUNCTION_CONFIRM_WINDOW)continue;

		int is_donor_test_ok=0;

		if(normally_arranged)
		{
			unsigned int small_pos_test_begin = left_virtualHead_abs_offset + (is_left_half_negative?real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW:(read_len - real_split_point_as_reversed)); 
			char * small_pos_read_begin = (is_left_half_negative?read_text:positive_read) + (is_left_half_negative?
						(real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW)           :
						(read_len - real_split_point_as_reversed)
  						);

			unsigned int large_pos_test_begin = right_virtualHead_abs_offset + (is_right_half_negative?real_split_point_as_reversed:(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));
			char * large_pos_read_begin = (is_right_half_negative?read_text:positive_read) + (is_right_half_negative?
						(real_split_point_as_reversed)     :
						(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));

			left_should_match = match_chro(small_pos_read_begin , value_index , small_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			right_should_match = match_chro(large_pos_read_begin , value_index , large_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			left_should_not_match = right_should_not_match = 0;
		//match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, left_virtualHead_abs_offset + real_split_point - JUNCTION_CONFIRM_WINDOW , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);

		}
		else
		{
			unsigned int small_pos_test_begin = left_virtualHead_abs_offset + (is_left_half_negative?real_split_point_as_reversed:(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW)); 
			char * small_pos_read_begin = (is_left_half_negative?read_text:positive_read) + (is_left_half_negative?
							(real_split_point_as_reversed):(read_len - real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW));

			unsigned int large_pos_test_begin = right_virtualHead_abs_offset + (is_right_half_negative?(real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW):(read_len - real_split_point_as_reversed));
			char * large_pos_read_begin = (is_right_half_negative?read_text:positive_read) + (is_right_half_negative?
							  (real_split_point_as_reversed - JUNCTION_CONFIRM_WINDOW):(read_len - real_split_point_as_reversed));

			left_should_match = match_chro(small_pos_read_begin , value_index , small_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			right_should_match = match_chro(large_pos_read_begin , value_index , large_pos_test_begin , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
			left_should_not_match = right_should_not_match = 0;

		}

		if(left_should_match + right_should_match  >= JUNCTION_CONFIRM_WINDOW*2 -1  &&
			left_should_not_match <= JUNCTION_CONFIRM_WINDOW -3 && right_should_not_match <= JUNCTION_CONFIRM_WINDOW -3)
		{
			int test_score = is_donor_test_ok*500+left_should_match + right_should_match - left_should_not_match - right_should_not_match;
			if(test_score > best_score)
			{
				selected_real_split_point = real_split_point_as_reversed;
				best_score = test_score;
			}
		}
	}

	if(best_score>0)
	{
		*final_split_point = selected_real_split_point;
		*is_donor_found = best_score>=500;
		*is_GT_AG_strand = selected_junction_strand;
		return best_score;
	}
	return 0;
}


int donor_score(global_context_t * global_context, thread_context_t * thread_context, unsigned int left_virtualHead_abs_offset, unsigned int right_virtualHead_abs_offset, int left_indel_offset, int right_indel_offset, int normally_arranged, int guess_start, int guess_end,  char * read_text, int read_len, int is_second_read, int * final_split_point, int * is_GT_AG_strand, int * is_donor_found)
{


	gene_value_index_t * value_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	int need_donor_test = global_context->config.is_rna_seq_reads;
	
	// guess_end is the index of the first UNWANTED BASE.
	int most_likely_point = (guess_start+guess_end)/2;
	
	// "split_point" is the first base NOT IN piece 1; it is also the first base IN piece 2. 
	int selected_real_split_point = -1, selected_junction_strand = -1;
	char donor_left[3], donor_right[3];
	
 
	int best_score = -111111;

	int real_split_point_i;
	int real_split_point_numbers = guess_end - guess_start;

	//printf("TESTDON: LR=%d; RR=%d\n", left_indel_offset, right_indel_offset);
	
	for(real_split_point_i = 0 ; real_split_point_i < real_split_point_numbers; real_split_point_i++)
	{
		int left_should_match, right_should_match;
		int left_should_not_match, right_should_not_match;
		int real_split_point = (real_split_point_i % 2)?-((real_split_point_i+1)/2):((1+real_split_point_i)/2);
		real_split_point += most_likely_point;
		int is_donor_test_ok = 0;

		if(real_split_point > read_len-JUNCTION_CONFIRM_WINDOW)continue;
		if(real_split_point < JUNCTION_CONFIRM_WINDOW)continue;

		assert(left_virtualHead_abs_offset<right_virtualHead_abs_offset);
		
		if(normally_arranged)
		{
			gvindex_get_string (donor_left, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, 2, 0);
			gvindex_get_string (donor_right, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset - 2, 2, 0);
		}
		else
		{
			gvindex_get_string (donor_left, value_index, right_virtualHead_abs_offset + real_split_point + left_indel_offset, 2, 0);
			gvindex_get_string (donor_right, value_index, left_virtualHead_abs_offset + real_split_point + right_indel_offset - 2, 2, 0);
		}
		is_donor_test_ok = is_donor_chars(donor_left) && is_donor_chars(donor_right) && paired_chars(donor_left, donor_right,0);


		donor_left[2]=0; donor_right[2]=0;
	//	printf("TESTDON: %s %s; OFFSET=%d; DON_OK=%d; NORMAL=%d; LEFT_OFF=%d; RIGHT_OFF=%d\n", donor_left, donor_right, real_split_point_i, is_donor_test_ok, normally_arranged, left_indel_offset, right_indel_offset);

		if(is_donor_test_ok || !need_donor_test)
		{
			if(normally_arranged)
			{
				left_should_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, left_virtualHead_abs_offset + real_split_point - JUNCTION_CONFIRM_WINDOW + left_indel_offset , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
				right_should_match = match_chro(read_text + real_split_point, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	

				left_should_not_match = match_chro(read_text + real_split_point, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
				right_should_not_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, right_virtualHead_abs_offset  + real_split_point + right_indel_offset - JUNCTION_CONFIRM_WINDOW, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
			}
			else
			{
				right_should_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, right_virtualHead_abs_offset + right_indel_offset + real_split_point - JUNCTION_CONFIRM_WINDOW , JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);
				left_should_match = match_chro(read_text + real_split_point, value_index, left_virtualHead_abs_offset + real_split_point + left_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	

				right_should_not_match = match_chro(read_text + real_split_point, value_index, right_virtualHead_abs_offset + real_split_point + right_indel_offset, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
				left_should_not_match = match_chro(read_text + real_split_point - JUNCTION_CONFIRM_WINDOW, value_index, left_virtualHead_abs_offset + left_indel_offset + real_split_point - JUNCTION_CONFIRM_WINDOW, JUNCTION_CONFIRM_WINDOW , 0, global_context -> config.space_type);	
	
			}

		//printf("!! TESTDON: M=%d,%d MM=%d,%d\n", left_should_match,+right_should_match,left_should_not_match,right_should_not_match);
			if(left_should_match +right_should_match >= 2*JUNCTION_CONFIRM_WINDOW-1 && 
				left_should_not_match <= JUNCTION_CONFIRM_WINDOW -5 && right_should_not_match <= JUNCTION_CONFIRM_WINDOW -5)
			{
				int test_score = is_donor_test_ok*3000+left_should_match + right_should_match - left_should_not_match - right_should_not_match;
				if(test_score > best_score)
				{
					selected_junction_strand = (donor_left[0]=='G' || donor_right[1]=='G');
					selected_real_split_point = real_split_point;	
					best_score = test_score;
				}
			}
		}
	}
	if(best_score>0)
	{
		*final_split_point = selected_real_split_point;
		*is_donor_found = best_score>=2900;
		*is_GT_AG_strand = selected_junction_strand;
		return best_score;
	}
	return 0;

}


void find_new_junctions(global_context_t * global_context, thread_context_t * thread_context, int pair_number, char * read_text, char * qual_text, int read_len, int is_second_read, int best_read_id)
{
	alignment_result_t * result =_global_retrieve_alignment_ptr(global_context, pair_number, is_second_read, best_read_id);
	subjunc_result_t * subjunc_result =_global_retrieve_subjunc_ptr(global_context, pair_number, is_second_read, best_read_id);


	if(read_len > EXON_LONG_READ_LENGTH)
		core_search_short_exons(global_context, thread_context,  read_text, qual_text, read_len, result -> selected_position, (subjunc_result -> minor_votes < 1)? result -> selected_position:subjunc_result -> minor_position, result -> confident_coverage_start, result -> confident_coverage_end);

	if(subjunc_result -> minor_votes < 1)return;
	if(result -> selected_votes < global_context->config.minimum_subread_for_first_read)return;

	if(global_context->config.do_big_margin_filtering_for_junctions)
	{
	//	if(2999302633 == result -> selected_position)
	//		printf("P0\n");
		if(is_ambiguous_voting(global_context, pair_number, is_second_read, result->selected_votes, result -> confident_coverage_start, result -> confident_coverage_end, read_len, (result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0))return;
	//	if(2999302633 == result -> selected_position)
	//		printf("P1\n");
	}

	/*
	if(2999302633 == result -> selected_position)
	{
		printf("MAIN_POS=%u; MINOR_POS=%u\n", result -> selected_position, subjunc_result -> minor_position);
		printf("SPLIT=%d\n",  subjunc_result->split_point);
	}*/

	unsigned int left_virtualHead_abs_offset = min(result -> selected_position, subjunc_result -> minor_position);
	unsigned int right_virtualHead_abs_offset = max(result -> selected_position, subjunc_result -> minor_position);

	int selected_real_split_point = subjunc_result->split_point;
	int is_GT_AG_donors = result->result_flags & 0x3;
	int is_donor_found = is_GT_AG_donors<3;
	int is_strand_jumped = (result->result_flags & CORE_IS_STRAND_JUMPED)?1:0;

	if(selected_real_split_point>0)
	{
		unsigned int left_edge_wanted, right_edge_wanted;

		if(is_strand_jumped)
		{
		
			// recover the "negative view" splicing point location
			int S = (result->result_flags & CORE_IS_NEGATIVE_STRAND) ? selected_real_split_point : (read_len - selected_real_split_point);
			int Sbar = read_len - S;

			int is_abnormal_as_reversed = (subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) == 1;
			if(!(result->result_flags & CORE_IS_NEGATIVE_STRAND)) is_abnormal_as_reversed = !is_abnormal_as_reversed;
			int is_small_half_negative = ((result->result_flags & CORE_IS_NEGATIVE_STRAND)?1:0) + (subjunc_result->minor_position < result->selected_position) ==1;

			if(is_abnormal_as_reversed && is_small_half_negative)
			{
				left_edge_wanted = left_virtualHead_abs_offset + S;
				right_edge_wanted = right_virtualHead_abs_offset + Sbar;
			}
			else if(is_abnormal_as_reversed && !is_small_half_negative)
			{
				left_edge_wanted = left_virtualHead_abs_offset + Sbar - 1;
				right_edge_wanted = right_virtualHead_abs_offset + S - 1;
			}
			else if(!is_abnormal_as_reversed && is_small_half_negative)
			{
				left_edge_wanted = left_virtualHead_abs_offset + S - 1;
				right_edge_wanted = right_virtualHead_abs_offset + Sbar - 1;
			}
			else // if(!is_abnormal_as_reversed && !is_small_half_negative)
			{
				left_edge_wanted = left_virtualHead_abs_offset + Sbar;
				right_edge_wanted = right_virtualHead_abs_offset + S;
			}
		}
		else
		{
			int selected_real_split_point_for_left = selected_real_split_point;
			int selected_real_split_point_for_right = selected_real_split_point;
			if((subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) == 1) //abnormally arranged halves
				selected_real_split_point_for_right --;
			else	// normally arranged halves
				selected_real_split_point_for_left --;


			
			int minor_indel_offset = (subjunc_result->double_indel_offset & 0xf);
			int major_indel_offset = (subjunc_result->double_indel_offset >> 4) & 0xf;
			if(major_indel_offset>=8)major_indel_offset=-(16-major_indel_offset);
			//assert(minor_indel_offset==0);
			//assert(major_indel_offset==0);

			left_edge_wanted = left_virtualHead_abs_offset + selected_real_split_point_for_left + ((result -> selected_position > subjunc_result -> minor_position)?minor_indel_offset: major_indel_offset);
			right_edge_wanted = right_virtualHead_abs_offset + selected_real_split_point_for_right;
		}
			

		//insert event
		HashTable * event_table = NULL;
		chromosome_event_t * event_space = NULL;
		if(thread_context)
		{
			event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
			event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
		}
		else
		{
			event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
			event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
		}

		// note that selected_real_split_point is the first UNWANTED base after left half.

		chromosome_event_t * found = NULL;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		int found_events = search_event(global_context, event_table, event_space, left_edge_wanted , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == right_edge_wanted)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found) found -> supporting_reads ++;
		else
		{
			int event_no;


			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;


			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = left_edge_wanted;
			new_event -> event_large_side = right_edge_wanted;


			if(is_donor_found &&(!is_strand_jumped) && right_edge_wanted - left_edge_wanted <= global_context -> config.maximum_intron_length
				&& (subjunc_result->minor_coverage_start > result->confident_coverage_start) + (subjunc_result -> minor_position >  result -> selected_position) !=1)
			{
				new_event -> is_negative_strand= !is_GT_AG_donors;
				new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

				new_event -> supporting_reads = 1;
				new_event -> indel_length = 0;
				
				put_new_event(event_table, new_event , event_no);

			}
			else
			{
				if(global_context -> config.do_fusion_detection)
				{
					new_event -> event_type = CHRO_EVENT_TYPE_FUSION;
					new_event -> is_strand_jumped = is_strand_jumped;

					new_event -> supporting_reads = 1;
					new_event -> indel_length = 0;
					
					put_new_event(event_table, new_event , event_no);
				}
			}
		}
	}
}

int write_fusion_final_results(global_context_t * global_context)
{
	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char fn2 [MAX_FILE_NAME_LENGTH];

	snprintf(fn2, MAX_FILE_NAME_LENGTH, "%s.fuse", global_context->config.output_prefix);
	FILE * ofp = f_subr_open(fn2, "wb");

	int xk1;
	//unsigned int all_junctions = 0;
	int no_sup_juncs = 0;
	int all_juncs = 0;

	for(xk1 = 0; xk1 < indel_context -> total_events ; xk1++)
	{ 
		char * chro_name_left,* chro_name_right;
		unsigned int chro_pos_left, chro_pos_right; 
		chromosome_event_t * event_body = indel_context -> event_space_dynamic +xk1;
		if(event_body -> event_type != CHRO_EVENT_TYPE_FUSION)
			continue;

		all_juncs++;
		if(event_body->final_counted_reads<1)
		{
			no_sup_juncs++;
			continue;
		}
		locate_gene_position( event_body -> event_small_side , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position( event_body -> event_large_side , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);

		chro_pos_left++;

		fprintf(ofp, "%s\t%u\t%s\t%u\t%c\t%d\n", chro_name_left, chro_pos_left, chro_name_right, chro_pos_right, event_body -> is_strand_jumped?'X':'=', event_body -> final_counted_reads);
	}

	fclose(ofp);
	return 0;
}
int write_junction_final_results(global_context_t * global_context)
{

	int no_sup_juncs = 0;

	indel_context_t * indel_context = (indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]; 
	char fn2 [MAX_FILE_NAME_LENGTH];

	snprintf(fn2, MAX_FILE_NAME_LENGTH, "%s.bed", global_context->config.output_prefix);
	FILE * ofp = f_subr_open(fn2, "wb");

	int xk1;
	unsigned int all_junctions = 0;

	for(xk1 = 0; xk1 < indel_context -> total_events ; xk1++)
	{ 
		char * chro_name_left,* chro_name_right;
		unsigned int chro_pos_left, chro_pos_right; 
		chromosome_event_t * event_body = indel_context -> event_space_dynamic +xk1;
		if(event_body -> event_type != CHRO_EVENT_TYPE_JUNCTION)
			continue;
		if(event_body->final_counted_reads<1)
		{
			no_sup_juncs++;
			continue;
		}

		locate_gene_position( event_body -> event_small_side , &global_context -> chromosome_table, &chro_name_left, &chro_pos_left);
		locate_gene_position( event_body -> event_large_side , &global_context -> chromosome_table, &chro_name_right, &chro_pos_right);

		chro_pos_left++;

		unsigned int feature_start = max(0, chro_pos_left - event_body -> junction_flanking_left );
		unsigned int feature_end = chro_pos_right + event_body -> junction_flanking_right;

		all_junctions ++;

		fprintf(ofp,"%s\t%u\t%u\tJUNC%08u\t%d\t%c\t%u\t%u\t%d,0,%d\t2\t%d,%d\t0,%u\n", chro_name_left, feature_start,  feature_end,
												all_junctions,  event_body -> final_counted_reads, event_body->is_negative_strand?'-':'+',
												feature_start,  feature_end, event_body->is_negative_strand?0:255, event_body->is_negative_strand?255:0,
												 event_body -> junction_flanking_left, event_body -> junction_flanking_right, feature_end-feature_start-event_body -> junction_flanking_right);
	
	}

	fclose(ofp);
	global_context -> all_junctions = all_junctions;
	//printf("Non-support juncs=%d;  Final juncs = %d\n", no_sup_juncs, all_junctions);
	return 0;
}



void get_chro_2base(char *buf, gene_value_index_t * index, unsigned int pos, int is_negative_strand)
{
	gvindex_get_string (buf, index, pos, 2, is_negative_strand);
}


int paired_chars_part(char * ch1, char * ch2, int is_reverse)
{
	if (c2eq(ch1, ch2, "GT", "AG") || c2eq(ch1, ch2, "CT", "AC"))
	{
		if (is_reverse) if (ceq(ch1, "AG") || ceq(ch1, "AC")) return 1;
		if (!is_reverse) if (ceq(ch1, "CT") || ceq(ch1, "GT")) return 1;
	}
	return 0;
}
#define is_donar_chars_part(cc) (((cc)[0]=='G' && (cc)[1]=='T') || \
			    ((cc)[0]=='A' && (cc)[1]=='G') || \
			    ((cc)[0]=='A' && (cc)[1]=='C') || \
			    ((cc)[0]=='C' && (cc)[1]=='T')) 


#define SHORT_EXON_MIN_LENGTH 18
#define EXON_EXTENDING_SCAN 0
#define SHORT_EXON_WINDOW 6 
#define SHORT_EXON_EXTEND 5000

void core_search_short_exons(global_context_t * global_context, thread_context_t * thread_context, char * read_text, char * qualityb0, int rl, unsigned int P1_Pos, unsigned int P2_Pos, short read_coverage_start, short read_coverage_end)
{
	char inb[1201], qualityb[1201];
	if ( (rl <= EXON_LONG_READ_LENGTH ) && (!EXON_EXTENDING_SCAN)) return;
	//return;
	gene_value_index_t * base_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;
	assert(base_index!=NULL);
	//insert event
	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}

	strcpy(inb, read_text);
	strcpy(qualityb, qualityb0);

	unsigned int pos_small=min(P1_Pos, P2_Pos), pos_big = max(P1_Pos, P2_Pos);

	int max_score , test_score;
	unsigned int best_j1_edge=0 , best_j2_edge=0;
	int need_to_test = 0;

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE HEAD  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

	if (read_coverage_start  > SHORT_EXON_MIN_LENGTH)
	{
		max_score = -1;

		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb , SHORT_EXON_MIN_LENGTH , global_context->config.phred_score_format); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(need_check2)
			if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb, base_index, pos_small, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type))
				need_check2 = 0; 


		if(need_check2)
		{

			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb, base_index, pos_small + delta_pos, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = SHORT_EXON_MIN_LENGTH ; test_donor_pos < read_coverage_start ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_small + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	max_score = -999;
	int max_is_GTAG = 0;

	if(need_to_test)
	{
		unsigned int test_end = pos_small - SHORT_EXON_EXTEND;
		if(SHORT_EXON_EXTEND > pos_small) test_end = 0;

		unsigned int new_pos = pos_small-16;
		while(1)
		{
			new_pos = match_chro_range(inb,  base_index, new_pos, 7 , new_pos - test_end , SEARCH_BACK);
			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// new_pos is the new head position of the read.
			int splice_point;
			for(splice_point = SHORT_EXON_MIN_LENGTH; splice_point < read_coverage_start ; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				get_chro_2base(cc, base_index, pos_small + splice_point -2, 0);
				if(is_donar_chars_part(cc))
				{
					// <<< EXON---|CC2---INTRON---CC|---EXON
					get_chro_2base(cc2, base_index, new_pos + splice_point, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc2 , cc, 0)) 
					{
						int matched_in_exon_old = match_chro(inb + splice_point, base_index, pos_small + splice_point , SHORT_EXON_WINDOW , 0, global_context->config.space_type);
						int matched_in_exon_new = match_chro(inb, base_index, new_pos , splice_point, 0, global_context->config.space_type);

						
						test_score = 1000000+ (matched_in_exon_new )*10000  + matched_in_exon_old * 1000 + new_pos - test_end;
						if(test_score <= max_score) continue;
						max_score = test_score + 39999 ;

						if(matched_in_exon_new < splice_point || matched_in_exon_old < SHORT_EXON_WINDOW ) 
							continue;

						max_is_GTAG = (cc2[0]=='G' || cc2[1]=='G');
						//printf("EX CC=%s\tCC2=%s\tis_GTAG=%d\n",cc,cc2,max_is_GTAG);
						best_j1_edge = new_pos + splice_point - 1;
						best_j2_edge = pos_small + splice_point;
					}
				}
			}
		}
	}


	if(best_j1_edge>0)
	{
		int event_no;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		chromosome_event_t * found = NULL;

		int found_events = search_event(global_context, event_table, event_space, best_j1_edge , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == best_j2_edge)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found) found -> supporting_reads ++;
		else
		{
			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = best_j1_edge;
			new_event -> event_large_side = best_j2_edge;
			assert(best_j1_edge<best_j2_edge);

			new_event -> is_negative_strand= !max_is_GTAG;
			new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

			new_event -> supporting_reads = 1;
			new_event -> indel_length = 0;

			put_new_event(event_table, new_event , event_no);
		}
		//printf("FOUND NEW JUNCTION HEAD: %u - %u\n", best_j1_edge, best_j2_edge);
	}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
// SCAN TO THE TAIL  /////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

	need_to_test = 0;
	max_score = -999;


	if (read_coverage_end< rl - SHORT_EXON_MIN_LENGTH)
	{
		int need_check2 = 1;
		if(qualityb[0])
		{
			float head_quality = read_quality_score(qualityb + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , global_context->config.phred_score_format); 
			if(head_quality < 6 )
				need_check2 = 0;
		}


		if(SHORT_EXON_MIN_LENGTH *0.6 < match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH , SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type))
			need_check2 = 0; 
		if(need_check2)
		{
			int delta_pos, is_indel = 0;
			for(delta_pos=-3; delta_pos <=3; delta_pos ++)
			{
				if(match_chro(inb + rl - SHORT_EXON_MIN_LENGTH, base_index, pos_big + rl - SHORT_EXON_MIN_LENGTH + delta_pos, SHORT_EXON_MIN_LENGTH , 0, global_context->config.space_type) >= SHORT_EXON_MIN_LENGTH*.7)
				{
					is_indel = 1;
					break;
				}
			}
			// The head of the read is incorrect. Do we need to search a long way?
			// See if there is a donor in the head area.
			int test_donor_pos;
			char cc[3];
			cc[2]=0;

			if(!is_indel)
				for(test_donor_pos = read_coverage_end  ; test_donor_pos < rl ; test_donor_pos ++)
				{
					get_chro_2base(cc, base_index, pos_big + test_donor_pos, 0);
					if(is_donar_chars_part(cc))
					{
						need_to_test = 1;
						break;
					}
				}
		}
	}

	best_j1_edge = 0;
	max_is_GTAG = 0;

	if(need_to_test)
	{
		unsigned int test_end = pos_big + SHORT_EXON_EXTEND;
		if(test_end > base_index -> length + base_index -> start_point) test_end = base_index -> length + base_index -> start_point;

		unsigned int new_pos = pos_big +rl - SHORT_EXON_MIN_LENGTH +16;

		while(1)
		{
			new_pos = match_chro_range(inb + rl - SHORT_EXON_MIN_LENGTH,  base_index, new_pos, 7 , test_end - new_pos , SEARCH_FRONT);
			if(new_pos==0xffffffff) break;
			// There is an exact match. See if the donor/receptors are matched.
			// (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point) is the new exon start.

			int splice_point;
			for(splice_point = read_coverage_end ; splice_point < rl -  SHORT_EXON_MIN_LENGTH; splice_point ++)
			{
				char cc[3];
				cc[2]=0;
				char cc2[3];
				cc2[2]=0;

				unsigned int new_pos_tail = (new_pos + SHORT_EXON_MIN_LENGTH -rl + splice_point);

				get_chro_2base(cc, base_index, pos_big + splice_point, 0);
				if(is_donar_chars_part(cc))
				{
					get_chro_2base(cc2, base_index, new_pos_tail -2, 0);
					if(is_donar_chars_part(cc2) && paired_chars_part(cc , cc2, 0)) 
					{
						int matched_in_exon_new = match_chro(inb + splice_point, base_index, new_pos_tail , rl - splice_point , 0, global_context->config.space_type);
						int matched_in_exon_old = match_chro(inb + splice_point - SHORT_EXON_WINDOW , base_index, pos_big + splice_point - SHORT_EXON_WINDOW , SHORT_EXON_WINDOW, 0, global_context->config.space_type);

						test_score = 1000000+ (matched_in_exon_new)*10000 + matched_in_exon_old * 1000  + test_end - new_pos;
						if(test_score <= max_score) continue;
						max_score = test_score + 39999;

						if(matched_in_exon_new < (rl - splice_point) || matched_in_exon_old < SHORT_EXON_WINDOW)
							continue;

						// EXON ---|CC---INTRON---CC2|--- EXON >>>
						max_is_GTAG = (cc[0]=='G'|| cc[1]=='G');
						best_j1_edge = pos_big + splice_point - 1;
						best_j2_edge = new_pos_tail;
					}
				}
			}

		}
	}


	if(best_j1_edge>0)
	{
		int event_no;
		chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
		chromosome_event_t * found = NULL;

		int found_events = search_event(global_context, event_table, event_space, best_j1_edge , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

		if(found_events)
		{
			int kx1; 
			for(kx1 = 0; kx1 < found_events ; kx1++)
			{
				if(search_return[kx1] -> event_large_side == best_j2_edge)
				{
					found = search_return[kx1];	
					break;
				}
			}
		}

		if(found) found -> supporting_reads ++;
		else
		{
			if(thread_context)
				event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
			else
				event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;


			event_space = reallocate_event_space(global_context, thread_context, event_no);

			chromosome_event_t * new_event = event_space+event_no; 
			memset(new_event,0,sizeof(chromosome_event_t));
			new_event -> event_small_side = best_j1_edge;
			new_event -> event_large_side = best_j2_edge;
			assert(best_j1_edge<best_j2_edge);

			new_event -> is_negative_strand= !max_is_GTAG;
			new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

			new_event -> supporting_reads = 1;
			new_event -> indel_length = 0;

			put_new_event(event_table, new_event , event_no);
			//printf("FOUND NEW JUNCTION TAIL: %u - %u\n", best_j1_edge, best_j2_edge);
		}
	}
}









int core_select_best_matching_halves_maxone(global_context_t * global_context, gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2, gehash_data_t max_pos, gene_vote_number_t max_votes, short max_start, short max_end, short max_mask, char * max_indel_recorder, int* best_select_max_votes, int rl)
{
	int best_splicing_point = -1, i,j;
	char * best_chro_name, is_reversed;
	unsigned int best_chro_pos;
	int selected_max_votes = -1;


	is_reversed = (max_mask & IS_NEGATIVE_STRAND)?1:0;
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			char * chro_name;
			char is_partner_reversed;
			unsigned int chro_pos;

			int overlapped_len, overlap_start, overlap_end;
			// All logical conditions

			//if( (vote->votes[i][j] < vote-> coverage_start[i][j]) < 12 && (vote-> coverage_end[i][j] > rl - 12 )) continue;

			is_partner_reversed = (vote->masks [i][j] & IS_NEGATIVE_STRAND) ? 1:0;
			overlap_start = max(max_start , vote->coverage_start[i][j]);
			overlap_end   = min(max_end , vote->coverage_end[i][j]);
			overlapped_len =overlap_end - overlap_start;

			int coverage_len = max_end - max_start + vote->coverage_end[i][j] - vote->coverage_start[i][j];
			if (overlapped_len >0)coverage_len -= overlapped_len;
			//SUBREADprintf("MAX: %d-%d   OTHER %d-%d    COV=%d   OVLP=%d\n", max_start, max_end, vote->coverage_start[i][j], vote->coverage_end[i][j], coverage_len, overlapped_len);



			if(overlapped_len >=14)
				continue;

			long long int dist = vote->pos[i][j];
			dist -= max_pos;

			//SUBREADprintf ("D=%lld\n", abs(dist));
			if (abs(dist)<6)
				continue;

			int support_r1 = 1; 
			int support_r2 = 1;

			if (max_votes < support_r1 || vote->votes[i][j]<support_r2)
				continue;

			// Same chromosome
			if ((vote->coverage_start[i][j] < max_start) + is_reversed == 1)
			{
				locate_gene_position(max_pos + read_len, &(global_context -> chromosome_table) , &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] , &(global_context -> chromosome_table), &chro_name, &chro_pos);
			}else
			{
				locate_gene_position(max_pos , &(global_context -> chromosome_table), &best_chro_name, &best_chro_pos);
				locate_gene_position(vote->pos[i][j] +read_len, &(global_context -> chromosome_table), &chro_name, &chro_pos);
			}

			if (chro_name != best_chro_name)	// The pointers can be compared because they can be the same.
				continue;

			int is_fusion = 0;

			if(is_reversed != is_partner_reversed) is_fusion = 1; 

			if( is_reversed && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] < max_start) != 1))is_fusion = 1;
			if((! is_reversed) && ((max_pos > vote->pos[i][j]) + (vote->coverage_start[i][j] > max_start) != 1)) is_fusion = 1;

			if(abs(dist) > 500000 || chro_name != best_chro_name) continue;

			int test_vote_value ;
			test_vote_value = 8888888 +  vote->votes[i][j]* 1000000 - abs(dist);
			if (hint_pos>=0)
			{
				long long int hint_dist = hint_pos;
				hint_dist -= vote->pos[i][j];
				if (abs (hint_dist) < 100000)
					test_vote_value += 100;
				if (abs (hint_dist) < 5000)
					test_vote_value += 100;
			}

			if (test_vote_value<selected_max_votes)continue;
			// Conditions of order of R3 and R5
			*half_marks &= ~IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] < max_start && (((max_pos < vote->pos[i][j]) && !is_reversed) || ((max_pos > vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;
			if (vote->coverage_start[i][j] >= max_end  &&  (((max_pos > vote->pos[i][j]) && !is_reversed) || ((max_pos < vote->pos[i][j]) && is_reversed) ) )
				*half_marks |= IS_REVERSED_HALVES;

			if (vote->coverage_start[i][j] < max_start)
			{
				(*half_marks) = (*half_marks) & ~IS_R1_CLOSE_TO_5;
			}
			else
			{
				(*half_marks) |= IS_R1_CLOSE_TO_5;
			}

			if(max_mask & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R1;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R1;

			if(vote->masks[i][j] & IS_NEGATIVE_STRAND)
				*half_marks = (*half_marks) |   IS_NEGATIVE_STRAND_R2;
			else
				*half_marks = (*half_marks) &  ~IS_NEGATIVE_STRAND_R2;
	

			
			best_splicing_point = ((vote->coverage_start[i][j] < max_start)? (vote->coverage_end[i][j]):(max_end)) + ((vote->coverage_start[i][j] < max_start)? (max_start):(vote->coverage_start[i][j]));


			best_splicing_point /=2;

			* best_pos1 = max_pos ;
			* best_pos2 = vote->pos[i][j] ;
			* best_vote1 = max_votes ;
			* best_vote2 = vote->votes[i][j] ;
			* read_coverage_start = min(vote->coverage_start[i][j] , max_start);
			* read_coverage_end = max(vote->coverage_end[i][j] , max_end);

			* read_coverage_start = max_start;
			* read_coverage_end = max_end;
			
			int k;
			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!max_indel_recorder[k+3])break;
			* indel_in_p1 = max_indel_recorder[k+2];

			for(k=0; k<MAX_INDEL_TOLERANCE ; k+=3)
				if(!vote->indel_recorder[i][j][k+3])break;
			* indel_in_p2 = vote->indel_recorder[i][j][k+2];


			* is_reversed_halves = is_reversed;

			if (test_vote_value >=100)
				*half_marks = (*half_marks) | IS_PAIRED_HINTED;
			else
				*half_marks = (*half_marks) & ~(IS_PAIRED_HINTED);

			if (is_fusion)
				*half_marks = (*half_marks)    | IS_FUSION;
			else
				*half_marks = (*half_marks) & ~( IS_FUSION);
	

			selected_max_votes = test_vote_value; 

		}
	*best_select_max_votes = selected_max_votes ;
	return best_splicing_point;
}



int core_select_best_matching_halves(global_context_t * global_context , gene_vote_t * vote, unsigned int * best_pos1, unsigned int * best_pos2, int * best_vote1, int * best_vote2, char * is_abnormal, short * half_marks, int * is_reversed_halves, float accept_rate, int read_len, long long int hint_pos, int tolerable_bases, short * read_coverage_start, short * read_coverage_end, char * indel_in_p1, char * indel_in_p2 , int * max_cover_start, int * max_cover_end, int rl, int repeated_pos_base, int is_negative, char * repeat_record, unsigned int index_valid_range)
{
	unsigned int tmp_best_pos1=0, tmp_best_pos2=0;
	int tmp_best_vote1=0, tmp_best_vote2=0, tmp_is_reversed_halves=0;
	char tmp_is_abnormal=0, tmp_indel_in_p1=0, tmp_indel_in_p2=0;
	short tmp_half_marks=0, tmp_read_coverage_start=0, tmp_read_coverage_end=0;
	int ret = 0, best_ret = 0;	

	int i,j;
	int test_select_votes=-1, best_select_votes = 1000000;
	//int max_minor = 0;

	/*
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor >= vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					max_minor = vote->votes[ii][jj];
				}

		}

	int encountered = 0;


	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(vote->votes[i][j] < vote->max_vote)continue;
			int ii,jj;
			for (ii=0; ii<GENE_VOTE_TABLE_SIZE;ii++)
				for(jj=0; jj< vote->items[ii]; jj++)
				{
					if(max_minor != vote->votes[ii][jj]) continue;
					if(ii==i && jj==j)continue;
					long long int dist =  vote->pos[ii][jj];
					dist =abs(dist - vote->pos[i][j]);
					if(dist > 500000)
						continue;
					encountered++;
				}

		}
	*/

	int repeated_pos = repeated_pos_base;
	int offset_shifting = (rl > 220)?4:0;
	//int encounter = 0;

	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			/*if((vote->votes[i][j] >=  vote->max_vote -1) && (vote->max_coverage_start >= vote-> coverage_start[i][j] - EXON_MAX_BIGMARGIN_OVERLAPPING ) &&  (vote->max_coverage_end <= vote-> coverage_end[i][j] + EXON_MAX_BIGMARGIN_OVERLAPPING))
				encounter++;*/
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] >=  vote->max_vote && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}
	for (i=0; i<GENE_VOTE_TABLE_SIZE; i++)
		for(j=0; j< vote->items[i]; j++)
		{
			if(repeated_pos_base>=0 && vote->pos[i][j]<=index_valid_range)
				if(vote->votes[i][j] ==  vote->max_vote -1 && repeated_pos < repeated_pos_base+12)
				{
					repeat_record[repeated_pos] = (vote-> coverage_start[i][j] >> offset_shifting);
					repeat_record[repeated_pos+1] = (vote-> coverage_end[i][j] >> offset_shifting);
					repeat_record[repeated_pos+2] = (is_negative?0x80:0) | (vote->votes[i][j]&0x7f);
					repeated_pos+=3;
				}
		}


	/*
	if(encounter>=2)
		return 0;
	*/

	ret = core_select_best_matching_halves_maxone(global_context, vote, &tmp_best_pos1, &tmp_best_pos2, &tmp_best_vote1, &tmp_best_vote2,  &tmp_is_abnormal,&tmp_half_marks, &tmp_is_reversed_halves, accept_rate, read_len, hint_pos,  tolerable_bases, &tmp_read_coverage_start, &tmp_read_coverage_end, &tmp_indel_in_p1, &tmp_indel_in_p2, vote -> max_position,  vote->max_vote, vote-> max_coverage_start, vote-> max_coverage_end,  vote-> max_mask, vote->max_indel_recorder, &test_select_votes, rl);
	test_select_votes += vote->max_vote*1000000;
			//SUBREADprintf("TSV=%d\n",test_select_votes);

	if(test_select_votes > best_select_votes)
	{
		best_select_votes = test_select_votes;
		*best_pos1 = tmp_best_pos1;
		*best_pos2 = tmp_best_pos2;
		*is_reversed_halves= tmp_is_reversed_halves;
		
		*best_vote1 = tmp_best_vote1;
		*best_vote2 = tmp_best_vote2;
		*is_abnormal = tmp_is_abnormal;
		*indel_in_p1 = tmp_indel_in_p1;
		*indel_in_p2 = tmp_indel_in_p2;
				
		*half_marks = tmp_half_marks;
		*read_coverage_start = tmp_read_coverage_start;
		*read_coverage_end = tmp_read_coverage_end;

		* max_cover_start = vote-> max_coverage_start;
		* max_cover_end = vote-> max_coverage_end;
		best_ret = ret;
	}		
	return best_ret;
}



#define EXON_DONOR_TEST_WINDOW 17


// pos1 must be small than pos2.
int core13_test_donor(char *read, int read_len, unsigned int pos1, unsigned int pos2, int guess_break_point, char negative_strand, int test_range, char is_soft_condition, int EXON_INDEL_TOLERANCE, int* real_break_point, gene_value_index_t * my_value_array_index, int indel_offset1, int indel_offset2, int is_reversed, int space_type, int * best_donor_score, int * is_GTAG)
{
	int bps_pos_x;
	int search_start = guess_break_point - test_range ;
	int search_end   = guess_break_point + test_range ;
	char h1_2ch[3], h2_2ch[3];

	h1_2ch[2] = h2_2ch[2]=0;
	search_start=max(10, search_start);
	search_end = min(read_len-10, search_end);
	int best_break = -1;
	int min_x = -9099;

	for (bps_pos_x = search_start; bps_pos_x < search_end ; bps_pos_x ++)
	{
		int paired_score = 0;
		get_chro_2base(h1_2ch, my_value_array_index, pos1 - indel_offset1+ bps_pos_x , is_reversed);
		get_chro_2base(h2_2ch, my_value_array_index, pos2 - 2 - indel_offset2 + bps_pos_x, is_reversed);


		//if(!is_reversed)
		//SUBREADprintf("C1=%s @%u, C2=%s @%u\n",h1_2ch, pos1 + bps_pos_x, h2_2ch,pos2 - 2 + indel_offset + bps_pos_x);
		if(h1_2ch[0]==h2_2ch[0] && h1_2ch[1]==h2_2ch[1]) continue;

		if(is_donar_chars_part(h1_2ch) && is_donar_chars_part(h2_2ch))
		{

			paired_score = paired_chars_part(h1_2ch, h2_2ch, is_reversed);

			if(paired_score)
			{
				int m1, m2, x1, x2;
				int break_point_half = is_reversed?(read_len - bps_pos_x):bps_pos_x;
				int first_exon_end,second_half_start;
				int donar_conf_len = 0;

				donar_conf_len = min(break_point_half , EXON_DONOR_TEST_WINDOW);
				donar_conf_len = min(read_len - break_point_half, donar_conf_len);
				//SUBREADprintf("DONOR_CONF_LEN=%d\n", donar_conf_len);

				if (is_reversed)
				{
					first_exon_end = pos2 + bps_pos_x - indel_offset2;
					second_half_start = pos1 + bps_pos_x- indel_offset1;

					m1 = match_chro(read + break_point_half - donar_conf_len , my_value_array_index, first_exon_end, donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start-donar_conf_len , donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end - donar_conf_len, donar_conf_len , is_reversed, space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len ,  my_value_array_index, second_half_start , donar_conf_len, is_reversed, space_type);
				}
				else
				{
					first_exon_end = pos1 + bps_pos_x - indel_offset1;
					second_half_start = pos2 + bps_pos_x - indel_offset2;

					m1 = match_chro(read + break_point_half - donar_conf_len, my_value_array_index, first_exon_end-donar_conf_len , donar_conf_len, is_reversed, space_type);
					m2 = match_chro(read + break_point_half , my_value_array_index, second_half_start, donar_conf_len, is_reversed, space_type);

					x1 = match_chro(read + break_point_half ,  my_value_array_index, first_exon_end, donar_conf_len , is_reversed,space_type);
					x2 = match_chro(read + break_point_half - donar_conf_len,  my_value_array_index, second_half_start - donar_conf_len, donar_conf_len , is_reversed,space_type);
				}

				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
				{
					SUBREADprintf("DONOR TEST STR=%s, %s ; pos=%d    %d %d ; M=%d %d ; X=%d %d\n", h1_2ch, h2_2ch, bps_pos_x, indel_offset1, indel_offset2, m1, m2, x1, x2);
				}
				#endif
	
				int threshold = 3;
				if (paired_score == 1)
					threshold = 3;

				#ifdef QUALITY_KILL
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len-1 )
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len- threshold )
				#else
				if (m1 >= donar_conf_len-1    && m2>=donar_conf_len -1)
					if(x1<donar_conf_len - threshold  && x2<donar_conf_len - threshold)
				#endif
					{
						int score =  3000-(x1 + x2) + (m1+ m2) ;
						if (min_x < score)
						{
							min_x = score;
							best_break = bps_pos_x;
							*is_GTAG = 1==((is_reversed) + (h1_2ch[0]=='G' || h1_2ch[1]=='G'));	//"GT" or "AG"
							//printf("FL CC=%s\tCC2=%s\tis_GTAG=%d\tREV=%d\n",h1_2ch,h2_2ch,*is_GTAG, is_reversed);
							*best_donor_score = score;
						}
					}
			}
		}
	}

	if (best_break>0)
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("SELECRED!!!_BREAKPOINT=%d, RAW POS=%u,%u, R=%s\n",  best_break, pos1 , pos2, read);
				#endif
		//SUBREADprintf ("FINAL BREAK: %d   ; REV = %d\n ", best_break, is_reversed);
		*real_break_point = best_break;
		return 1;
	}
	else
	{
				#ifdef TEST_TARGET
				if(memcmp(read, TEST_TARGET, 15)==0)
					SUBREADprintf("KILLED!!!_BREAKPOINT=%d, R=%s\n",  best_break+ pos1, read);
				#endif
	}
	return 0;
}






#define EXON_LARGE_WINDOW 60
#define ACCEPTED_SUPPORT_RATE 0.3

void core_fragile_junction_voting(global_context_t * global_context, thread_context_t * thread_context, char * read, char * qual, unsigned int full_rl, int negative_strand, int color_space, unsigned int low_border, unsigned int high_border, gene_vote_t *vote_p1)
{
	int windows = full_rl / EXON_LARGE_WINDOW +1;
	float overlap = (1.0*windows * EXON_LARGE_WINDOW - full_rl) / (windows-1);

	int ww;
	int window_cursor = 0;

	HashTable * event_table = NULL;
	chromosome_event_t * event_space = NULL;
	if(thread_context)
	{
		event_table = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}
	else
	{
		event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_entry_table; 
		event_space = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) -> event_space_dynamic;
	}



	for(ww=0; ww<windows;ww++)
	{
		window_cursor = (int)(ww * EXON_LARGE_WINDOW - ww * overlap);
		int read_len = EXON_LARGE_WINDOW;
		if(ww == windows-1)
			read_len = full_rl -window_cursor;

		float subread_step = 3.00001;
		int i;
		int subread_no;
		char * InBuff;
		InBuff = read + window_cursor;
		char tmp_char = InBuff[read_len];
		InBuff[read_len] = 0;
		
		init_gene_vote(vote_p1);
		for(subread_no=0; ; subread_no++)
		{
			int subread_offset1 = (int)(subread_step * (subread_no+1));
			subread_offset1 -= subread_offset1%GENE_SLIDING_STEP;
			subread_offset1 += GENE_SLIDING_STEP-1;

			for(i=0; i<GENE_SLIDING_STEP ; i++)
			{
				int subread_offset = (int)(subread_step * subread_no); 
				subread_offset -= subread_offset%GENE_SLIDING_STEP -i;

				char * subread_string = InBuff + subread_offset;
				gehash_key_t subread_integer = genekey2int(subread_string, color_space);

				gehash_go_q(global_context->current_index, subread_integer , subread_offset, read_len,negative_strand, vote_p1, 1, 1, 21.9, 24, 5, subread_no,  low_border, high_border - read_len);
			}
			if(subread_offset1 >= read_len -16)
				break;
		}


		if(1)
		{
			finalise_vote(vote_p1);
			select_best_vote(vote_p1);
			//print_votes(vote_p1, global_context -> config.index_prefix);
			unsigned int best_pos1=0;
			unsigned int best_pos2=0;
			int best_vote1=0;
			int best_vote2=0;
			char is_abnormal=0;
			short half_marks=0;
			int is_reversed_halves=0, max_cover_start=0, max_cover_end=0;
			char indel_in_p1=0, indel_in_p2=0;
			short read_coverage_start =0, read_coverage_end=0;
			gene_value_index_t * base_index = thread_context?thread_context->current_value_index:global_context->current_value_index ;

			int splice_point = core_select_best_matching_halves(global_context, vote_p1, &best_pos1, &best_pos2, &best_vote1, &best_vote2, &is_abnormal ,&half_marks, &is_reversed_halves, ACCEPTED_SUPPORT_RATE, read_len, -1,  0, &read_coverage_start, &read_coverage_end, &indel_in_p1, &indel_in_p2, &max_cover_start, &max_cover_end, read_len, -1 , 0, NULL , 0xffffffff);

			//printf("SP=%d;  BV=%d;  BV2=%d\n", splice_point, best_vote1, best_vote2);
			if (splice_point>0 && best_vote1 >= 1 && best_vote2>=1)
			{
				int test_real_break_point = -1, test_donor_score=-1;
				int is_GTAG = 0;
				int is_accepted = core13_test_donor(InBuff, read_len, min(best_pos1, best_pos2), max(best_pos1,best_pos2), splice_point, negative_strand, read_len/4, 0, 5, &test_real_break_point, base_index, 0, 0, negative_strand, color_space, &test_donor_score, &is_GTAG);

				if (is_accepted ){
					unsigned int pos_small = min(test_real_break_point+ best_pos1,  test_real_break_point+ best_pos2) - 1;
					unsigned int pos_big = max(test_real_break_point+ best_pos1,  test_real_break_point+ best_pos2);

					int event_no;
					chromosome_event_t * search_return [MAX_EVENT_ENTRIES_PER_SITE];
					chromosome_event_t * found = NULL;

					int found_events = search_event(global_context, event_table, event_space, pos_small , EVENT_SEARCH_BY_SMALL_SIDE, CHRO_EVENT_TYPE_JUNCTION|CHRO_EVENT_TYPE_FUSION, search_return);

					if(found_events)
					{
						int kx1; 
						for(kx1 = 0; kx1 < found_events ; kx1++)
						{
							if(search_return[kx1] -> event_large_side == pos_big)
							{
								found = search_return[kx1];	
								break;
							}
						}
					}

					if(found) found -> supporting_reads ++;
					else
					{
						if(thread_context)
							event_no = ((indel_thread_context_t *)thread_context -> module_thread_contexts[MODULE_INDEL_ID]) -> total_events ++;
						else
							event_no = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID]) ->  total_events ++;

						event_space = reallocate_event_space(global_context, thread_context, event_no);

						chromosome_event_t * new_event = event_space+event_no; 
						memset(new_event,0,sizeof(chromosome_event_t));
						new_event -> event_small_side = pos_small;
						new_event -> event_large_side = pos_big;

						new_event -> is_negative_strand= !is_GTAG;
						new_event -> event_type = CHRO_EVENT_TYPE_JUNCTION;

						new_event -> supporting_reads = 1;
						new_event -> indel_length = 0;

						put_new_event(event_table, new_event , event_no);
						//printf("ADD JUNCTION BY FRAGILE, %d-%d\n", pos_small, pos_big);
					}


				}

			}
		}
		InBuff[read_len] = tmp_char;
	}
}


