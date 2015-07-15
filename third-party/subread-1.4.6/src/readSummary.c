/***************************************************************

   The Subread and Rsubread software packages are free
   software packages:
 
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
  
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>


#ifndef MAKE_STANDALONE
  #include <R.h>
#endif

#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <getopt.h>
#include "subread.h"
#include "interval_merge.h"
#include "core.h"
#include "gene-algorithms.h"
#include "sambam-file.h"
#include "input-files.h"
#include "hashtable.h"
#include "HelperFunctions.h"

/********************************************************************/
/********************************************************************/
/********************************************************************/
//  NEW FUNCTION FOR MULTI-THREADING
/********************************************************************/
/********************************************************************/
/********************************************************************/
#define FEATURE_NAME_LENGTH  256 
#define CHROMOSOME_NAME_LENGTH 256 
#define MAX_LINE_LENGTH 3000
#define FILE_TYPE_RSUBREAD 10
#define FILE_TYPE_GTF 100

#define ALLOW_ALL_MULTI_MAPPING 1
#define ALLOW_PRIMARY_MAPPING 2

#define MAX_HIT_NUMBER 3000

typedef struct
{
	unsigned int feature_name_pos;
	unsigned int start;
	unsigned int end;
	unsigned int sorted_order;

	unsigned short chro_name_pos_delta;
	char is_negative_strand;
} fc_feature_info_t;

typedef struct
{
	unsigned long long assigned_reads;
	unsigned long long unassigned_ambiguous;
	unsigned long long unassigned_multimapping;
	unsigned long long unassigned_nofeatures;
	unsigned long long unassigned_unmapped;
	unsigned long long unassigned_mappingquality;
	unsigned long long unassigned_fragmentlength;
	unsigned long long unassigned_chimericreads;
	unsigned long long unassigned_secondary;
	unsigned long long unassigned_nonjunction;
	unsigned long long unassigned_duplicate;
} fc_read_counters;

typedef struct
{
	unsigned short thread_id;
	char * line_buffer1;
	char * line_buffer2;
	unsigned long long int nreads_mapped_to_exon;
	unsigned long long int all_reads;
	//unsigned short current_read_length1;
	//unsigned short current_read_length2;
	unsigned int * count_table;
	unsigned int chunk_read_ptr;
	pthread_t thread_object;

	char * input_buffer;
	unsigned int input_buffer_remainder;
	unsigned int input_buffer_write_ptr;	
	pthread_spinlock_t input_buffer_lock;


	short hits_total_length1[MAX_HIT_NUMBER];
	short hits_total_length2[MAX_HIT_NUMBER];
	long hits_indices1 [MAX_HIT_NUMBER];
	long hits_indices2 [MAX_HIT_NUMBER];
	long decision_table_ids [MAX_HIT_NUMBER];
	unsigned char decision_table_votes [MAX_HIT_NUMBER];
	long decision_table_exon_ids [MAX_HIT_NUMBER];
	long uniq_gene_exonid_table [MAX_HIT_NUMBER];
	long uniq_gene_table [MAX_HIT_NUMBER];

	char * chro_name_buff;
	z_stream * strm_buffer;

	fc_read_counters read_counters;

	SamBam_Alignment aln_buffer;
} fc_thread_thread_context_t;

#define REVERSE_TABLE_BUCKET_LENGTH 131072
#define REDUCE_TO_5_PRIME_END 5
#define REDUCE_TO_3_PRIME_END 3

typedef struct
{
	unsigned int chro_number;
	unsigned int chro_features;
	unsigned int chro_feature_table_start;
	unsigned int chro_block_table_start;
	unsigned int chro_block_table_end;
	unsigned int chro_possible_length;

	unsigned short chro_reverse_table_current_size;
	unsigned int * reverse_table_start_index;
	//unsigned int * reverse_table_end_index;
} fc_chromosome_index_info;

typedef struct
{
	int is_gene_level;
	int is_paired_end_data;
	int is_multi_overlap_allowed;
	int is_strand_checked;
	int is_both_end_required;
	int is_chimertc_disallowed;
	int is_PE_distance_checked;
	int is_multi_mapping_allowed;
	int is_input_file_resort_needed;
	int is_SAM_file;
	int is_read_details_out;
	int is_unpaired_warning_shown;
	int is_stake_warning_shown;
	int is_split_alignments_only;
	int is_duplicate_ignored;
	int reduce_5_3_ends_to_one;
	int isCVersion;

	int min_mapping_quality_score;
	int min_paired_end_distance;
	int max_paired_end_distance;
	int feature_block_size;
	int read_length;
	int line_length;
	int longest_chro_name;
	int five_end_extension;
	int three_end_extension;
	int overlap_length_required;

	unsigned long long int all_reads;

	unsigned short thread_number;
	fc_thread_thread_context_t * thread_contexts;
	int is_all_finished;
	unsigned int input_buffer_max_size;
	SamBam_Reference_Info * sambam_chro_table;

	char * debug_command;
	char * unistr_buffer_space;
	unsigned int unistr_buffer_size;
	unsigned int unistr_buffer_used;

	HashTable * gene_name_table;	// gene_name -> gene_number
	HashTable * annot_chro_name_alias_table;	// name in annotation file -> alias name
	char alias_file_name[300];
	char input_file_name[300];
	char raw_input_file_name[300];
	char output_file_name[300];
	unsigned char ** gene_name_array;	// gene_internal_number -> gene_name 

	HashTable * exontable_chro_table;	// gene_name -> fc_chromosome_index_info structure (contains chro_number, feature_number, block_start, block_end, etc) 
	int exontable_nchrs;
	int exontable_exons;
	int * exontable_geneid;
	char * exontable_strand;
	char ** exontable_chr;
	long * exontable_start;
	long * exontable_stop;
	char feature_name_column[100];
	char gene_id_column[100];

	long * exontable_block_end_index;
	long * exontable_block_max_end;
	long * exontable_block_min_start;

	char ** exontable_anno_chrs;
	char * exontable_anno_chr_2ch;
	long * exontable_anno_chr_heads;

	FILE * SAM_output_fp;
	double start_time;

	char * cmd_rebuilt;

	char   redo;

	fc_read_counters read_counters;
	
} fc_thread_global_context_t;

unsigned int tick_time = 1000;


unsigned int unistr_cpy(fc_thread_global_context_t * global_context, char * str, int strl)
{
	unsigned int ret;
	if(global_context->unistr_buffer_used + strl >= global_context->unistr_buffer_size-1)
	{
		if( global_context->unistr_buffer_size < 3435973835u) // 4G / 5 * 4 - 5
		{
			global_context -> unistr_buffer_size = global_context->unistr_buffer_size /4 *5;
			global_context -> unistr_buffer_space = realloc(global_context -> unistr_buffer_space, global_context->unistr_buffer_size);
		}
		else
		{
			SUBREADprintf("Error: exceed memory limit (4GB) for storing annotation data.\n");
			return 0xffffffffu;
		}
	}

	strcpy(global_context -> unistr_buffer_space + global_context->unistr_buffer_used, str);
	ret = global_context->unistr_buffer_used;

	global_context->unistr_buffer_used += strl +1;

	return ret;
}

void print_FC_configuration(fc_thread_global_context_t * global_context, char * annot, char * sam, char * out, int is_sam, int is_GTF, int *n_input_files, int isReadSummaryReport)
{
	char * tmp_ptr1 = NULL , * next_fn, *sam_used = malloc(strlen(sam)+1), sam_ntxt[30],bam_ntxt[30], next_ntxt[50];
	int nfiles=1, nBAMfiles = 0, nNonExistFiles = 0;

	strcpy(sam_used, sam);

	SUBREADputs("");
	print_subread_logo();
	SUBREADputs("");
	print_in_box(80,1,1,"featureCounts setting");
	print_in_box(80,0,0,"");
	
	nfiles = 0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		nfiles++;

		int file_probe = is_certainly_bam_file(next_fn, NULL);
		if(file_probe==-1) nNonExistFiles++;
		if(file_probe == 1) nBAMfiles++;		
	}

	sam_ntxt[0]=0;
	bam_ntxt[0]=0;
	next_ntxt[0]=0;

	if(nNonExistFiles)
		sprintf(next_ntxt, "%d unknown file%s", nNonExistFiles, nNonExistFiles>1?"s":"");
	if(nBAMfiles)
		sprintf(bam_ntxt, "%d BAM file%s  ", nBAMfiles, nBAMfiles>1?"s":"");
	if(nfiles-nNonExistFiles-nBAMfiles)
		sprintf(sam_ntxt, "%d SAM file%s  ", nfiles-nNonExistFiles-nBAMfiles , (nfiles-nNonExistFiles-nBAMfiles)>1?"s":"");


	strcpy(sam_used, sam);

	print_in_box(80,0,0,"            Input files : %s%s%s", sam_ntxt, bam_ntxt, next_ntxt);
	nfiles=0;

	while(1)
	{
		next_fn = strtok_r(nfiles==0?sam_used:NULL, ";", &tmp_ptr1);
		if(next_fn == NULL || strlen(next_fn)<1) break;
		int is_first_read_PE = 0 , file_probe = is_certainly_bam_file(next_fn, &is_first_read_PE);

		char file_chr = 'S';
		if(file_probe == -1) file_chr = '?';
		else if(is_first_read_PE == 1) file_chr = 'P';
		//file_chr = 'o';

		print_in_box(94,0,0,"                          %c[32m%c%c[36m %s%c[0m",CHAR_ESC, file_chr,CHAR_ESC, next_fn,CHAR_ESC);
		nfiles++;
	}

	(*n_input_files) = nfiles;
	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"            Output file : %s", out);
	print_in_box(80,0,0,"            Annotations : %s (%s)", annot, is_GTF?"GTF":"SAF");
	if(isReadSummaryReport)
		print_in_box(80,0,0,"     Assignment details : <input_file>.featureCounts");

	if(global_context -> alias_file_name[0])
		print_in_box(80,0,0,"  Chromosome alias file : %s", global_context -> alias_file_name);

	print_in_box(80,0,0,"");
	print_in_box(80,0,0,"                Threads : %d", global_context->thread_number);
	print_in_box(80,0,0,"                  Level : %s level", global_context->is_gene_level?"meta-feature":"feature");
	print_in_box(80,0,0,"             Paired-end : %s", global_context->is_paired_end_data?"yes":"no");
	print_in_box(80,0,0,"        Strand specific : %s", global_context->is_strand_checked?(global_context->is_strand_checked==1?"yes":"inversed"):"no");
	char * multi_mapping_allow_mode = "not counted";
	if(global_context->is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING)
		multi_mapping_allow_mode = "primary only";
	else if(global_context->is_multi_mapping_allowed == ALLOW_ALL_MULTI_MAPPING)
		multi_mapping_allow_mode = "counted";
	print_in_box(80,0,0,"     Multimapping reads : %s", multi_mapping_allow_mode);
	print_in_box(80,0,0,"Multi-overlapping reads : %s", global_context->is_multi_overlap_allowed?"counted":"not counted");
	if(global_context -> is_split_alignments_only)
		print_in_box(80,0,0,"       Split alignments : required");
	if(global_context -> overlap_length_required !=1)
		print_in_box(80,0,0,"      Overlapping bases : %d", global_context -> overlap_length_required);
	if(global_context -> five_end_extension || global_context -> three_end_extension)
		print_in_box(80,0,0,"        Read extensions : %d on 5' and %d on 3' ends", global_context -> five_end_extension , global_context -> three_end_extension);
	if(global_context -> reduce_5_3_ends_to_one)
		print_in_box(80,0,0,"      Read reduction to : %d' end" , global_context -> reduce_5_3_ends_to_one == REDUCE_TO_5_PRIME_END ?5:3);
	if(global_context -> is_duplicate_ignored)
		print_in_box(80,0,0,"       Duplicated Reads : ignored");

	if(global_context->is_paired_end_data)
	{
		print_in_box(80,0,0,"");
		print_in_box(80,0,0,"         Chimeric reads : %s", global_context->is_chimertc_disallowed?"not counted":"counted");
		print_in_box(80,0,0,"       Both ends mapped : %s", global_context->is_both_end_required?"required":"not required");

		if(global_context->is_PE_distance_checked)
			print_in_box(80,0,0,"        Fragment length : %d - %d", global_context -> min_paired_end_distance, global_context -> max_paired_end_distance);
	}

	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	print_in_box(80,1,1,"Running");
	print_in_box(80,0,0,"");
	if(global_context->annot_chro_name_alias_table)
		print_in_box(80,0,0,"%ld chromosome name aliases are loaded.", global_context -> annot_chro_name_alias_table ->numOfElements);

	free(sam_used);
}

void print_FC_results(fc_thread_global_context_t * global_context)
{
	print_in_box(89,0,1,"%c[36mRead assignment finished.%c[0m", CHAR_ESC, CHAR_ESC);
	print_in_box(80,0,0,"");
	print_in_box(80,2,1,"http://subread.sourceforge.net/");
	SUBREADputs("");
	return;


	if(0){
		print_in_box(80,1,1,"Summary");
		print_in_box(80,0,0,"");
		if(global_context->is_paired_end_data)
			print_in_box(80,0,0,"        All fragments : %llu", global_context -> all_reads);
		else
			print_in_box(80,0,0,"            All reads : %llu", global_context -> all_reads);

		if(global_context->is_gene_level)
			print_in_box(80,0,0,"        Meta-features : %lu", global_context -> gene_name_table -> numOfElements);
		else
			print_in_box(80,0,0,"             Features : %u", global_context -> exontable_exons);

		if(global_context->is_paired_end_data)
			print_in_box(80,0,0,"   Assigned fragments : %llu", global_context -> read_counters.assigned_reads);
		else
			print_in_box(80,0,0,"       Assigned reads : %llu", global_context -> read_counters.assigned_reads);

		print_in_box(80,0,0,"            Time cost : %.3f minutes", (miltime() - global_context -> start_time)/60);
		print_in_box(80,0,0,"");
		print_in_box(80,2,1,"http://subread.sourceforge.net/");
	}
	SUBREADputs("");
}

int fc_strcmp(const void * s1, const void * s2)
{
	return strcmp((char*)s1, (char*)s2);
}


int is_comment_line(const char * l, int file_type, unsigned int lineno)
{
	int tabs = 0, xk1 = 0;
	if(l[0]=='#') return 1;

	if(isalpha(l[0]) && file_type == FILE_TYPE_RSUBREAD)
	{
		char target_chr[16];
		memcpy(target_chr, l, 16);
		for(xk1=0; xk1<16; xk1++)
			target_chr[xk1] = tolower(target_chr[xk1]);

		if(memcmp(target_chr, "geneid\tchr\tstart",16)==0) return 1;
	}

	xk1=0;
	while(l[xk1]) tabs += (l[xk1++] == '\t');

	return tabs < ((file_type == FILE_TYPE_GTF)?8:4);
}

// This function loads annotations from the file.
// It returns the number of featres loaded, or -1 if something is wrong. 
// Memory will be allowcated in this function. The pointer is saved in *loaded_features.
// The invoker must release the memory itself.

int load_feature_info(fc_thread_global_context_t *global_context, const char * annotation_file, int file_type, fc_feature_info_t ** loaded_features)
{
	unsigned int features = 0, xk1 = 0, lineno=0;
	char * file_line = malloc(MAX_LINE_LENGTH+1);
	FILE * fp = f_subr_open(annotation_file,"r"); 
	int is_GFF_warned = 0;
	if(!fp) return -1;

	HashTable * chro_name_table = HashTableCreate(1603);
	HashTableSetHashFunction(chro_name_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(chro_name_table, fc_strcmp_chro);
	global_context -> longest_chro_name = 0;

	// first scan: get the chromosome size, etc
	while(1)
	{
		char * fgets_ret = fgets(file_line, MAX_LINE_LENGTH, fp);
		char * token_temp, *chro_name;
		fc_chromosome_index_info * chro_stab;
		unsigned int feature_pos = 0;
		if(!fgets_ret) break;

		lineno++;
		if(is_comment_line(file_line, file_type, lineno-1))continue;
		if(file_type == FILE_TYPE_GTF)
		{
			chro_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp); // lib_name (not needed)
			char * feature_type = strtok_r(NULL,"\t", &token_temp);
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				strtok_r(NULL,"\t", &token_temp); // feature_start
				feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end
				features++;
			}
			else chro_name = NULL;
		}
		else
		{
			strtok_r(file_line,"\t", &token_temp);
			chro_name = strtok_r(NULL,"\t",&token_temp);
			strtok_r(NULL,"\t",&token_temp);	// feature_start
			feature_pos = atoi(strtok_r(NULL,"\t", &token_temp));// feature_end

			features++;
		}

		if(chro_name)
		{
			if(strlen(chro_name)>=CHROMOSOME_NAME_LENGTH) 
				chro_name[CHROMOSOME_NAME_LENGTH-1]=0;
			chro_stab = HashTableGet(chro_name_table, chro_name);

			if(chro_stab)
			{
				chro_stab -> chro_possible_length = max(chro_stab -> chro_possible_length , feature_pos+1);
			}else
			{
				char * tmp_chro_name = malloc(CHROMOSOME_NAME_LENGTH);
				term_strncpy(tmp_chro_name, chro_name, CHROMOSOME_NAME_LENGTH);
				chro_stab = calloc(sizeof(fc_chromosome_index_info),1);
				chro_stab -> chro_number = chro_name_table->numOfElements;
				chro_stab -> chro_possible_length = feature_pos+1;
				chro_stab -> reverse_table_start_index = NULL;
				HashTablePut(chro_name_table, tmp_chro_name, chro_stab);
			}

			chro_stab -> chro_features ++;
		}
	}

	fseek(fp,0,SEEK_SET);

	fc_feature_info_t * ret_features = malloc(sizeof(fc_feature_info_t) * features);

	lineno = 0;
	while(xk1 < features)
	{
		int is_gene_id_found = 0;
		fgets(file_line, MAX_LINE_LENGTH, fp);
		lineno++;
		char * token_temp;
		if(is_comment_line(file_line, file_type, lineno-1))continue;

		if(file_type == FILE_TYPE_RSUBREAD)
		{
			char * feature_name = strtok_r(file_line,"\t",&token_temp);
			int feature_name_len = strlen(feature_name);
			if(feature_name_len > FEATURE_NAME_LENGTH) feature_name[FEATURE_NAME_LENGTH -1 ] = 0;
			ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name, feature_name_len);

			char * seq_name = strtok_r(NULL,"\t", &token_temp);
			int chro_name_len = strlen(seq_name);
			if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
			unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
			global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

			ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;
			ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
			if(ret_features[xk1].start<0 || ret_features[xk1].start>0x7fffffff)
			{
				ret_features[xk1].start = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}

			ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 
			if(ret_features[xk1].end<0 || ret_features[xk1].end>0x7fffffff)
			{
				ret_features[xk1].end = 0;
				print_in_box(80,0,0,"WARNING the %d-th line has a negative chro coordinate.", lineno);
			}




			char * strand_str = strtok_r(NULL,"\t", &token_temp); 
			if(strand_str == NULL)
				ret_features[xk1].is_negative_strand = 0;
			else
				ret_features[xk1].is_negative_strand = ('-' ==strand_str[0]);
			ret_features[xk1].sorted_order = xk1;

			int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
			
			fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
			if(!chro_stab -> reverse_table_start_index)
			{
				chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
			}
			chro_stab -> reverse_table_start_index[bin_location]++;

			is_gene_id_found = 1;
			xk1++;
		}
		else if(file_type == FILE_TYPE_GTF)
		{
			char feature_name_tmp[FEATURE_NAME_LENGTH];
			sprintf(feature_name_tmp, "LINE_%07u", xk1 + 1);
			char * seq_name = strtok_r(file_line,"\t",&token_temp);
			strtok_r(NULL,"\t", &token_temp);// source
			char * feature_type = strtok_r(NULL,"\t", &token_temp);// feature_type
			if(strcmp(feature_type, global_context -> feature_name_column)==0)
			{
				ret_features[xk1].start = atoi(strtok_r(NULL,"\t", &token_temp));// start 
				ret_features[xk1].end = atoi(strtok_r(NULL,"\t", &token_temp));//end 

				if(ret_features[xk1].start < 1 || ret_features[xk1].end<1 ||  ret_features[xk1].start > 0x7fffffff ||  ret_features[xk1].end > 0x7fffffff || ret_features[xk1].start > ret_features[xk1].end)
					SUBREADprintf("\nWarning: the feature on the %d-th line has zero coordinate or zero lengths\n\n", lineno);


				strtok_r(NULL,"\t", &token_temp);// score 
				ret_features[xk1].is_negative_strand = ('-' == (strtok_r(NULL,"\t", &token_temp)[0]));//strand 
				ret_features[xk1].sorted_order = xk1;
				strtok_r(NULL,"\t",&token_temp);	// "frame"
				char * extra_attrs = strtok_r(NULL,"\t",&token_temp);	// name_1 "val1"; name_2 "val2"; ... 
				if(extra_attrs && (strlen(extra_attrs)>2))
				{
					int attr_val_len = GTF_extra_column_value(extra_attrs , global_context -> gene_id_column , feature_name_tmp, FEATURE_NAME_LENGTH);
					if(attr_val_len>0) is_gene_id_found=1;
			//		printf("V=%s\tR=%d\n", extra_attrs , attr_val_len);
				}

				if(is_gene_id_found)
				{
				}
				else
				{
					if(!is_GFF_warned)
					{
						int ext_att_len = strlen(extra_attrs);
						if(extra_attrs[ext_att_len-1] == '\n') extra_attrs[ext_att_len-1] =0;
						SUBREADprintf("\nWarning: failed to find the gene identifier attribute in the 9th column of the provided GTF file.\nThe specified gene identifier attribute is '%s' \nThe attributes included in your GTF annotation are '%s' \n\n",  global_context -> gene_id_column, extra_attrs);
					}
					is_GFF_warned++;
				}

				int feature_name_len = strlen(feature_name_tmp);
				if(feature_name_len > FEATURE_NAME_LENGTH) feature_name_tmp[FEATURE_NAME_LENGTH -1 ] = 0;
				ret_features[xk1].feature_name_pos = unistr_cpy(global_context, (char *)feature_name_tmp, feature_name_len);

				int chro_name_len = strlen(seq_name);
				if(chro_name_len > CHROMOSOME_NAME_LENGTH) seq_name[CHROMOSOME_NAME_LENGTH -1 ] = 0;
				unsigned int chro_name_pos = unistr_cpy(global_context, (char *)seq_name, chro_name_len);
				global_context -> longest_chro_name = max(chro_name_len, global_context -> longest_chro_name);

				ret_features[xk1].chro_name_pos_delta = chro_name_pos - ret_features[xk1].feature_name_pos;

				int bin_location = ret_features[xk1].start / REVERSE_TABLE_BUCKET_LENGTH;
				fc_chromosome_index_info * chro_stab = HashTableGet(chro_name_table, seq_name);
				if(!chro_stab -> reverse_table_start_index)
				{
					chro_stab -> reverse_table_start_index = malloc(sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
					memset(chro_stab -> reverse_table_start_index, 0 , sizeof(int) *( chro_stab->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2));
				}
				chro_stab -> reverse_table_start_index[bin_location]++;

				xk1++;
			}
		}
	}
	fclose(fp);
	free(file_line);

	(*loaded_features) = ret_features;
	global_context -> exontable_nchrs = (int)chro_name_table-> numOfElements;
	global_context -> exontable_chro_table = chro_name_table;

	print_in_box(80,0,0,"   Features : %d\n", features);
	if(features < 1)
	{
		print_in_box(80,0,0,"WARNING no features were loaded in format %s.", file_type == FILE_TYPE_GTF?"GTF":"SAF");
		print_in_box(80,0,0,"        annotation format can be specified using '-F'.");
	}
	return features;
}

int find_or_insert_gene_name(fc_thread_global_context_t * global_context, unsigned char * feature_name)
{
	HashTable * genetable = global_context -> gene_name_table;

	long long int gene_number = HashTableGet(genetable, feature_name) - NULL;
	if(gene_number>0)
		return gene_number-1;
	else
	{
		gene_number = genetable -> numOfElements; 
		HashTablePut(genetable, feature_name, NULL+gene_number+1);
		global_context -> gene_name_array[gene_number] = feature_name;
			// real memory space of feature_name is in the "loaded_features" data structure.
			// now we only save its pointer.

		return gene_number;
	}
}

void register_reverse_table(int block_no, long this_block_min_start, long this_block_max_end, fc_chromosome_index_info * chro_inf)
{

	unsigned int reversed_bucket_start = this_block_min_start /  REVERSE_TABLE_BUCKET_LENGTH;
	unsigned int reversed_bucket_end = this_block_max_end / REVERSE_TABLE_BUCKET_LENGTH;
	assert(this_block_min_start <= this_block_max_end);
	assert(reversed_bucket_end < chro_inf -> chro_possible_length);
	int x1;
	for(x1 = reversed_bucket_start; x1 <= reversed_bucket_end; x1++)
	{
		chro_inf->reverse_table_start_index[x1] = min(chro_inf->reverse_table_start_index[x1], block_no);
		//chro_inf->reverse_table_end_index[x1] = max(chro_inf->reverse_table_end_index[x1], block_no+1);
	}

}

void feature_merge(void * arrv, int start, int items, int items2)
{

	void ** arr = (void **) arrv;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	int total_items = items+items2;
	long * tmp_start = malloc(sizeof(long) * total_items);
	long * tmp_end = malloc(sizeof(long) * total_items);
	unsigned char * tmp_strand = malloc(sizeof(char) * total_items);
	int * tmp_entyrez = malloc(sizeof(int) * total_items);
	fc_feature_info_t ** tmp_info_ptr = malloc(sizeof(fc_feature_info_t*) * total_items);

	int read_1_ptr = start;
	int read_2_ptr = start+items;
	int write_ptr;

	for(write_ptr=0; write_ptr<total_items; write_ptr++)
	{
		if((read_1_ptr >= start+items)||(read_2_ptr < start+total_items && ret_start[read_1_ptr] >= ret_start[read_2_ptr]))
		{
			tmp_start[write_ptr] = ret_start[read_2_ptr];
			tmp_end[write_ptr] = ret_end[read_2_ptr];
			tmp_strand[write_ptr] = ret_strand[read_2_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_2_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_2_ptr];
			read_2_ptr++;
		}
		else
		{
			tmp_start[write_ptr] = ret_start[read_1_ptr];
			tmp_end[write_ptr] = ret_end[read_1_ptr];
			tmp_strand[write_ptr] = ret_strand[read_1_ptr];
			tmp_entyrez[write_ptr] = ret_entyrez[read_1_ptr];
			tmp_info_ptr[write_ptr] = old_info_ptr[read_1_ptr];
			read_1_ptr++;
		}
	}

	memcpy(ret_start+ start, tmp_start, sizeof(long) * total_items);
	memcpy(ret_end+ start, tmp_end, sizeof(long) * total_items);
	memcpy(ret_strand+ start, tmp_strand, sizeof(char) * total_items);
	memcpy(ret_entyrez+ start, tmp_entyrez, sizeof(int) * total_items);
	memcpy(old_info_ptr+ start, tmp_info_ptr, sizeof(fc_feature_info_t*) * total_items);

	free(tmp_start);
	free(tmp_end);
	free(tmp_strand);
	free(tmp_entyrez);
	free(tmp_info_ptr);
}


int feature_sort_compare(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long * ret_start = (long *)arr[0];
	long ll = ret_start[l];
	long rl = ret_start[r];

	if(ll==rl) return 0;
	else if(ll>rl) return 1;
	else return -1;
}

void feature_sort_exchange(void * arrv, int l, int r)
{
	void ** arr = (void **) arrv;
	long tmp;
	fc_feature_info_t * tmpptr;

	long * ret_start = (long *) arr[0];
	long * ret_end = (long *) arr[1];
	unsigned char * ret_strand = (unsigned char *) arr[2];
	int * ret_entyrez = (int *) arr[3];
	fc_feature_info_t ** old_info_ptr = (fc_feature_info_t **) arr[4];

	
	tmp = ret_start[r];
	ret_start[r]=ret_start[l];
	ret_start[l]=tmp;

	tmp = ret_end[r];
	ret_end[r]=ret_end[l];
	ret_end[l]=tmp;

	tmp = ret_strand[r];
	ret_strand[r]=ret_strand[l];
	ret_strand[l]=tmp;

	tmp = ret_entyrez[r];
	ret_entyrez[r]=ret_entyrez[l];
	ret_entyrez[l]=tmp;

	tmpptr = old_info_ptr[r];
	old_info_ptr[r]=old_info_ptr[l];
	old_info_ptr[l]=tmpptr;

}



void sort_feature_info(fc_thread_global_context_t * global_context, unsigned int features, fc_feature_info_t * loaded_features, char *** sorted_chr_names, int ** sorted_entrezid, long ** sorted_start, long ** sorted_end, unsigned char ** sorted_strand, char ** anno_chr_2ch, char *** anno_chrs, long ** anno_chr_head, long ** block_end_index, long ** block_min_start_pos, long ** block_max_end_pos)
{
	unsigned int chro_pnt;
	unsigned int xk1,xk2;
	int * ret_entrez = malloc(sizeof(int) * features);
	long * ret_start = malloc(sizeof(long) * features);
	long * ret_end = malloc(sizeof(long) * features);
	int current_block_buffer_size = 2000;

	long * ret_block_end_index = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_min_start = malloc(sizeof(long) * current_block_buffer_size);
	long * ret_block_max_end = malloc(sizeof(long) * current_block_buffer_size);
	unsigned char * ret_strand = malloc(features);
	char ** ret_char_name = malloc(sizeof(void *) * features);
	fc_feature_info_t ** old_info_ptr = malloc(sizeof(void *) * features);
	(*anno_chrs) = malloc(sizeof(void *) * global_context -> exontable_nchrs);
	(*anno_chr_head) = malloc(sizeof(long) * global_context -> exontable_nchrs);
	(*anno_chr_2ch) = malloc(sizeof(char) * global_context -> exontable_nchrs*2); 
	unsigned int * chro_feature_ptr = calloc(sizeof(int) * global_context -> exontable_nchrs,1);
	fc_chromosome_index_info ** tmp_chro_info_ptrs = malloc(global_context -> exontable_nchrs * sizeof(fc_chromosome_index_info *));

	global_context -> gene_name_array = malloc(sizeof(char *) * features);	// there should be much less identical names.
	global_context -> gene_name_table = HashTableCreate(5000);
	HashTableSetHashFunction(global_context -> gene_name_table, HashTableStringHashFunction);
	HashTableSetKeyComparisonFunction(global_context -> gene_name_table, fc_strcmp);

	// init start positions of each chromosome block.
	if(1)
	{
		KeyValuePair * cursor;
		int bucket;
		unsigned int sum_ptr = 0;
		for(bucket=0; bucket < global_context -> exontable_chro_table  -> numOfBuckets; bucket++)
		{
			cursor = global_context -> exontable_chro_table -> bucketArray[bucket];
			while (1)
			{
				if (!cursor) break;
				fc_chromosome_index_info * tmp_chro_inf = cursor -> value;
				cursor = cursor->next;
				//tmp_chro_inf -> reverse_table_end_index = calloc(sizeof(int), tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
				chro_feature_ptr [tmp_chro_inf -> chro_number] = tmp_chro_inf -> chro_features;
				tmp_chro_info_ptrs[tmp_chro_inf -> chro_number] = tmp_chro_inf;
			}
		}

		for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
		{
			unsigned int tmpv = chro_feature_ptr[xk1];
			chro_feature_ptr[xk1] = sum_ptr;
			tmp_chro_info_ptrs[xk1] -> chro_feature_table_start = sum_ptr;
		//		printf("SII=%u  +  %u\n", sum_ptr, tmpv);
			sum_ptr += tmpv;
		}

	}
	int current_block_id = 0, sort_i = 0;

	(*sorted_chr_names) = ret_char_name;
	(*sorted_entrezid) = ret_entrez;
	(*sorted_start) = ret_start;
	(*sorted_end) = ret_end;
	(*sorted_strand) = ret_strand;
	int curr_chro_number = 0;

	for(chro_pnt=0; chro_pnt < features; chro_pnt++)
	{
		char * this_chro_name = global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos + loaded_features[chro_pnt].chro_name_pos_delta;
		fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table , this_chro_name);
		assert(this_chro_info);
		unsigned int this_chro_number = this_chro_info -> chro_number;
		unsigned int this_chro_table_ptr = chro_feature_ptr[this_chro_number];

		ret_char_name[this_chro_table_ptr] = this_chro_name;// (char *)loaded_features[chro_pnt].chro;
		ret_entrez[this_chro_table_ptr] = find_or_insert_gene_name(global_context, (unsigned char *)(global_context -> unistr_buffer_space + loaded_features[chro_pnt].feature_name_pos));
		ret_start[this_chro_table_ptr] = loaded_features[chro_pnt].start;
		ret_end[this_chro_table_ptr] = loaded_features[chro_pnt].end;
		ret_strand[this_chro_table_ptr] = loaded_features[chro_pnt].is_negative_strand;
		old_info_ptr[this_chro_table_ptr] = &loaded_features[chro_pnt];

		chro_feature_ptr[this_chro_number]++;
	}

	for(xk1 = 0; xk1 < global_context -> exontable_nchrs; xk1++)
	{
		fc_chromosome_index_info * tmp_chro_inf = tmp_chro_info_ptrs[xk1];
		int bins_in_chr = ( tmp_chro_inf->chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH +2);
		short * features_per_block_bins = malloc(sizeof(short)*bins_in_chr);
		for(xk2=0; xk2<bins_in_chr; xk2++)
		{
			features_per_block_bins[xk2] = max(1,min(1000,(int)(0.9999999+sqrt(tmp_chro_inf -> reverse_table_start_index[xk2]))));
			//printf("CHR%d : SQR[%d]=%d (%d)\n",  tmp_chro_inf -> chro_number,xk2, features_per_block_bins[xk2], tmp_chro_inf -> reverse_table_start_index[xk2] );
		}

		memset(tmp_chro_inf -> reverse_table_start_index, 0xff, sizeof(int) *bins_in_chr);

		tmp_chro_inf -> chro_block_table_start = current_block_id; 
		unsigned int this_block_items = 0;
		long this_block_min_start = 0x7fffffff, this_block_max_end = 0;
		unsigned int this_chro_tab_end =  tmp_chro_inf -> chro_features + tmp_chro_inf -> chro_feature_table_start;

		void * in_array[5];
		in_array[0] = ret_start + tmp_chro_inf -> chro_feature_table_start; 
		in_array[1] = ret_end + tmp_chro_inf -> chro_feature_table_start; 
		in_array[2] = ret_strand + tmp_chro_inf -> chro_feature_table_start; 
		in_array[3] = ret_entrez + tmp_chro_inf -> chro_feature_table_start; 
		in_array[4] = old_info_ptr + tmp_chro_inf -> chro_feature_table_start; 

		merge_sort(in_array, this_chro_tab_end - tmp_chro_inf -> chro_feature_table_start, feature_sort_compare, feature_sort_exchange, feature_merge);

		for(sort_i = tmp_chro_inf -> chro_feature_table_start; sort_i< this_chro_tab_end ; sort_i++)
		{
			// NOW THE FEATURES (ret_start, ret_end, ret_strand, ret_entrez, old_info_ptr) ARE ALL SORTED!
			//printf("NT=%lu\tCHRO=%d\n", ret_start[sort_i], tmp_chro_inf->chro_number);
			old_info_ptr[sort_i]->sorted_order = sort_i;

			int feature_bin_location = ret_start[sort_i] / REVERSE_TABLE_BUCKET_LENGTH;
			int block_bin_location = this_block_min_start / REVERSE_TABLE_BUCKET_LENGTH;

			if(this_block_items && (this_block_items > features_per_block_bins[block_bin_location] || feature_bin_location != block_bin_location))//global_context -> feature_block_size)
			{

				if(current_block_id >= current_block_buffer_size - 1)
				{
					current_block_buffer_size *= 1.3;
					ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
					ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
					ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
				}


				ret_block_end_index[current_block_id] = sort_i;	// FIRST UNWANTED ID
				ret_block_min_start[current_block_id] = this_block_min_start;
				ret_block_max_end[current_block_id] = this_block_max_end;
				register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
				//printf("B=%d; ST=%ld, END=%ld, ITM=%d\n", current_block_id, this_block_min_start, this_block_max_end, this_block_items);
				current_block_id++;
				this_block_max_end = 0;
				this_block_items = 0;
				this_block_min_start = 0x7fffffff;
			}

			this_block_max_end = max(this_block_max_end, ret_end[sort_i]);
			this_block_min_start = min(this_block_min_start, ret_start[sort_i]);
			this_block_items ++;
		
		}
		if(this_block_items)
		{
			if(current_block_id >= current_block_buffer_size)
			{
				current_block_buffer_size *= 1.3;
				ret_block_min_start = realloc(ret_block_min_start, sizeof(long)*current_block_buffer_size);
				ret_block_max_end = realloc(ret_block_max_end, sizeof(long)*current_block_buffer_size);
				ret_block_end_index = realloc(ret_block_end_index, sizeof(long)*current_block_buffer_size);
			}

			ret_block_end_index[current_block_id] = this_chro_tab_end;	// FIRST UNWANTED ID
			ret_block_min_start[current_block_id] = this_block_min_start;
			ret_block_max_end[current_block_id] = this_block_max_end;
			register_reverse_table(current_block_id, this_block_min_start, this_block_max_end, tmp_chro_inf);
			current_block_id++;
		}

		(*anno_chr_head) [curr_chro_number] = current_block_id; 
		tmp_chro_inf -> chro_block_table_end = current_block_id; 
		free(features_per_block_bins);
	}

	(*block_end_index) = ret_block_end_index;
	(*block_min_start_pos) = ret_block_min_start;
	(*block_max_end_pos) = ret_block_max_end;

	//print_in_box(80, 0,0,"The %u features are sorted.\n", sort_i);
	free(old_info_ptr);
	free(tmp_chro_info_ptrs);
	free(chro_feature_ptr);
}

//#define MAX_HIT_NUMBER 1800
//
int strcmp_slash(char * s1, char * s2)
{
	char nch;
	while(0!=(nch = *(s1++))){
		if(nch == '/') break;
		if(nch != (*s2)) return 1;
		s2++;
	}
	return nch != *s2;
}


void process_line_buffer(fc_thread_global_context_t * global_context, fc_thread_thread_context_t * thread_context)
{

	char * read_chr, *read_1_chr = NULL, *tmp_tok_ptr= NULL, *CIGAR_str , *read_name = NULL, *read_name1 = NULL;
	long read_pos, fragment_length = 0, read_1_pos = 0;
	unsigned int search_start = 0, search_end;
	int nhits1 = 0, nhits2 = 0, alignment_masks, search_block_id, search_item_id;
	long * hits_indices1 = thread_context -> hits_indices1, * hits_indices2 = thread_context -> hits_indices2;
	short * hits_total_length1 = thread_context -> hits_total_length1 ,  * hits_total_length2 = thread_context -> hits_total_length2;

	int is_second_read;
	int skipped_for_exonic = 0;
	int first_read_quality_score = 0;

	thread_context->all_reads++;
	//if(thread_context->all_reads>1000000) printf("TA=%llu\n%s\n",thread_context->all_reads, thread_context -> line_buffer1);

	for(is_second_read = 0 ; is_second_read < 2; is_second_read++)
	{
		if(is_second_read && !global_context -> is_paired_end_data) break;

		char * line = is_second_read? thread_context -> line_buffer2:thread_context -> line_buffer1;
	
		//printf("LINE_BUF=%s\n",line);

		read_name = strtok_r(line,"\t", &tmp_tok_ptr);	// read name
		if(!read_name)return;

		if(is_second_read)
		{
			if(read_name)
			{
				int x1;
				for(x1=0; read_name[x1]; x1++)
				{
					if(read_name[x1]=='/')
					{
						read_name[x1]=0;
						read_name1[x1]=0;
						break;
					}
				}
				//printf("R1=%s; R2=%s\n",read_name,read_name1 );
				if(strcmp_slash(read_name,read_name1)!=0)
				{
					//printf("WARN:%d [%d]\n", global_context->is_unpaired_warning_shown, thread_context -> thread_id);
					if(!global_context->is_unpaired_warning_shown)
					{
				//	printf("RV:%s,%s\n", read_name, read_name1);
						global_context->is_unpaired_warning_shown=1;
						global_context->redo = 1;
						print_in_box(80,0,0,"   Found reads that are not properly paired.");
						print_in_box(80,0,0,"   (missing mate or the mate is not the next read)");
					}
				}
			}
		}
		else
				read_name1 = read_name;

		char * mask_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if((!mask_str) || !isdigit(mask_str[0])) return;

		alignment_masks = atoi(mask_str);

		if(is_second_read == 0)
		{
			//skip the read if unmapped (its mate will be skipped as well if paired-end)
			if( ((!global_context -> is_paired_end_data) &&  (alignment_masks & SAM_FLAG_UNMAPPED) ) ||
			    ((alignment_masks & SAM_FLAG_UNMAPPED)   &&  (alignment_masks & SAM_FLAG_MATE_UNMATCHED) && global_context -> is_paired_end_data) ||
			    (((alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED)) && global_context -> is_paired_end_data && global_context -> is_both_end_required)
			  ){
				thread_context->read_counters.unassigned_unmapped ++;

				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Unmapped\t*\t*\n", read_name);
				return;	// do nothing if a read is unmapped, or the first read in a pair of reads is unmapped.
			}
		}



		read_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!read_chr) return;
		char * read_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!read_pos_str) return;

		read_pos = atoi(read_pos_str);
		if(read_pos < 1 && read_pos_str[0]!='0') return;

		char * mapping_qual_str = strtok_r(NULL,"\t", &tmp_tok_ptr);

		CIGAR_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
		if(!CIGAR_str)
			continue;

		if(global_context -> min_mapping_quality_score>0)
		{
			int mapping_qual =atoi(mapping_qual_str);

			//printf("SECOND=%d; FIRST=%d; THIS=%d; Q=%d\n", is_second_read, first_read_quality_score, mapping_qual, );
			if(( mapping_qual < global_context -> min_mapping_quality_score  && ! global_context -> is_paired_end_data)||( is_second_read  && max( first_read_quality_score, mapping_qual ) < global_context -> min_mapping_quality_score))
			{
				thread_context->read_counters.unassigned_mappingquality ++;

				if(global_context -> SAM_output_fp)
				{
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MappingQuality\t*\tMapping_Quality=%d,%d\n", read_name, first_read_quality_score, mapping_qual);
				}
				return;
			}
			if(is_second_read==0 && global_context -> is_paired_end_data)
			{
				first_read_quality_score = mapping_qual;
			}
		}


		long mate_pos = 0;
		char * mate_chr = NULL;

		if(is_second_read)
		{
			mate_chr = strtok_r(NULL,"\t", &tmp_tok_ptr);// mate_chr
			if(mate_chr[0]=='=') mate_chr = read_chr;
			char * mate_pos_str = strtok_r(NULL,"\t", &tmp_tok_ptr);	// mate_pos
			mate_pos = atol(mate_pos_str);

		}

		if(is_second_read == 0 && global_context -> is_paired_end_data && 
	   	  (global_context -> is_PE_distance_checked || global_context -> is_chimertc_disallowed)
		  )
		{
			int is_half_mapped = (alignment_masks & SAM_FLAG_UNMAPPED) || (alignment_masks & SAM_FLAG_MATE_UNMATCHED);

			if(!is_half_mapped)
			{
				char * mate_chrx = strtok_r(NULL,"\t", &tmp_tok_ptr); //get chr which the mate read is mapped to
				if(!mate_chrx) return;
				strtok_r(NULL,"\t", &tmp_tok_ptr);
				if(!tmp_tok_ptr) return;
				char * frag_len_str = strtok_r(NULL,"\t", &tmp_tok_ptr);
				if(!tmp_tok_ptr) return;

				fragment_length = abs(atoi(frag_len_str)); //get the fragment length

				int is_first_read_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
				int is_second_read_negative_strand = (alignment_masks & SAM_FLAG_MATE_REVERSE_STRAND_MATCHED)?1:0; 

				if(mate_chrx[0]=='=' && is_first_read_negative_strand!=is_second_read_negative_strand)
				{
					if(global_context -> is_PE_distance_checked && ((fragment_length > global_context -> max_paired_end_distance) || (fragment_length < global_context -> min_paired_end_distance)))
					{
						thread_context->read_counters.unassigned_fragmentlength ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_FragmentLength\t*\tLength=%ld\n", read_name, fragment_length);
						return;
					}
				}
				else
				{
					if(global_context -> is_chimertc_disallowed)
					{
						thread_context->read_counters.unassigned_chimericreads ++;

						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Chimera\t*\t*\n", read_name);
						return;
					}
				}
			}
		}

		if(!tmp_tok_ptr) return;


		// This filter has to be put here because the 0x400 FLAG is not about mapping but about sequencing.
		// A unmapped read with 0x400 FLAG should be able to kill the mapped mate which may have no 0x400 FLAG. 
		if(global_context -> is_duplicate_ignored)
		{
			if(alignment_masks & SAM_FLAG_DUPLICATE)
			{
				thread_context->read_counters.unassigned_duplicate ++;
				if(global_context -> SAM_output_fp)
					fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Duplicate\t*\t*\n", read_name);

				return;
			}

		}

		if(SAM_FLAG_UNMAPPED & alignment_masks) continue;

		char * NH_pos = strstr(tmp_tok_ptr,"\tNH:i:");
		if(NH_pos)
		{
			if(NH_pos[6]>'1' || isdigit(NH_pos[7]))
			{

				if(is_second_read && read_1_chr)
				{
					if((strcmp(read_1_chr, mate_chr)!=0 || mate_pos!=read_1_pos) && read_1_chr[0] != '*'  && mate_chr[0]!='*')
					{
				//	printf("RV:%s,%s   %d,%d\n", read_1_chr, mate_chr, mate_pos, read_1_pos);
						if(!global_context->is_unpaired_warning_shown)
						{
							global_context->is_unpaired_warning_shown=1;
							global_context->redo = 1;
							print_in_box(80,0,0,"   Found reads that are not properly paired.");
							print_in_box(80,0,0,"   (missing mate or mate is not the next read)");
						}
					}
				}
				else
				{
					read_1_chr = read_chr;
					read_1_pos = read_pos;
				}


				if(global_context -> is_multi_mapping_allowed == 0)
				{
					// now it is a NH>1 read!
					// not allow multimapping -> discard!
					thread_context->read_counters.unassigned_multimapping ++;

					if(global_context -> SAM_output_fp)
						fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_MultiMapping\t*\t*\n", read_name);
					return;
				}
			}
		}

		// if a pair of reads have one secondary, the entire fragment is seen as secondary.
		if((alignment_masks & SAM_FLAG_SECONDARY_MAPPING) && (global_context -> is_multi_mapping_allowed == ALLOW_PRIMARY_MAPPING))
		{
			thread_context->read_counters.unassigned_secondary ++;

			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Secondary\t*\t*\n", read_name);
			return;
		}

		int is_this_negative_strand = (alignment_masks & SAM_FLAG_REVERSE_STRAND_MATCHED)?1:0; 
		int is_second_read_in_pair = alignment_masks & SAM_FLAG_SECOND_READ_IN_PAIR;
		int is_fragment_negative_strand = is_second_read_in_pair?(!is_this_negative_strand):is_this_negative_strand;

		fc_chromosome_index_info * this_chro_info = HashTableGet(global_context -> exontable_chro_table, read_chr);
		if(this_chro_info == NULL)
		{
			if(global_context -> annot_chro_name_alias_table)
			{
				char * anno_chro_name = HashTableGet( global_context -> annot_chro_name_alias_table , read_chr);
				if(anno_chro_name)
					this_chro_info = HashTableGet(global_context -> exontable_chro_table, anno_chro_name);
			}
			if(this_chro_info == NULL && memcmp(read_chr, "chr", 3)==0)
			{
				this_chro_info = HashTableGet(global_context -> exontable_chro_table, read_chr+3);
			}
			//printf("NL=%s, CI=%p\n", (read_chr), this_chro_info);
			if(this_chro_info == NULL && strlen(read_chr)<=2)
			{
				strcpy(thread_context -> chro_name_buff, "chr");
				strcpy(thread_context -> chro_name_buff+3, read_chr);
				this_chro_info = HashTableGet(global_context -> exontable_chro_table, thread_context -> chro_name_buff);
			}
		}

		if(this_chro_info)
		{
			int nhits = 0;

			int cigar_section_id, cigar_sections;
			unsigned int Starting_Points[6];
			unsigned short Section_Lengths[6];
			long * hits_indices = (is_second_read?hits_indices2:hits_indices1);
			short * hits_total_length = (is_second_read?hits_total_length2:hits_total_length1);

			cigar_sections = RSubread_parse_CIGAR_string(CIGAR_str, Starting_Points, Section_Lengths);


			if(cigar_sections>1 || !global_context->is_split_alignments_only) 
			{

			//#warning "=================== COMMENT THESE 2 LINES ================================"
			//for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			//	SUBREADprintf("BCCC: %llu , sec[%d] %d ~ %d ; secs=%d ; flags=%d ; second=%d\n", read_pos, cigar_section_id , Starting_Points[cigar_section_id], Section_Lengths[cigar_section_id], cigar_sections, alignment_masks, is_second_read);
				if(global_context -> reduce_5_3_ends_to_one)
				{
					if((REDUCE_TO_5_PRIME_END == global_context -> reduce_5_3_ends_to_one) + is_this_negative_strand == 1) // reduce to 5' end (small coordinate if positive strand / large coordinate if negative strand)
					{
						Section_Lengths[0]=1;
					}
					else
					{
						Starting_Points[0] = Starting_Points[cigar_sections-1] + Section_Lengths[cigar_sections-1] - 1;
						Section_Lengths[0]=1;
					}

					cigar_sections = 1;
				}

				// Extending the reads to the 3' and 5' ends. (from the read point of view) 
				if(global_context -> five_end_extension)
				{
					if(is_this_negative_strand){
						Section_Lengths [cigar_sections - 1] += global_context -> five_end_extension;
					}else{
						//SUBREADprintf("5-end extension: %d [%d]\n", Starting_Points[0], Section_Lengths[0]);
						if( read_pos > global_context -> five_end_extension)
						{
							Section_Lengths [0] += global_context -> five_end_extension;
							read_pos  -= global_context -> five_end_extension;
						}
						else
						{
							Section_Lengths [0] += read_pos-1;
							read_pos = 1;
						}
					}
				}

				if(global_context -> three_end_extension)
				{

					if(is_this_negative_strand){
						if( read_pos > global_context -> three_end_extension)
						{
							Section_Lengths [0] += global_context -> three_end_extension;
							read_pos -= global_context -> three_end_extension;
						}
						else
						{
							Section_Lengths [0] += read_pos - 1;
							read_pos = 1;
						}
					}
					else	Section_Lengths [cigar_sections - 1] += global_context -> three_end_extension;

				}

			//#warning "=================== COMMENT THESE 2 LINES ================================"
			//for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
			//	SUBREADprintf("ACCC: %llu , sec[%d] %u ~ %d ; secs=%d\n", read_pos, cigar_section_id, Starting_Points[cigar_section_id], Section_Lengths[cigar_section_id], cigar_sections);
				if( global_context -> debug_command [0] == 'D')
					SUBREADprintf("\n\nRead name = %s ; Second read = %d\n", read_name, is_second_read);

				for(cigar_section_id = 0; cigar_section_id<cigar_sections; cigar_section_id++)
				{
					long section_begin_pos = read_pos + Starting_Points[cigar_section_id];
					long section_end_pos = Section_Lengths[cigar_section_id] + section_begin_pos - 1;

					
					if(global_context -> debug_command [0] == 'D')
						SUBREADprintf("Section [%d]: %d ~ %d\n", cigar_section_id , Starting_Points[cigar_section_id], Starting_Points[cigar_section_id]+Section_Lengths[cigar_section_id]);

					int start_reverse_table_index = section_begin_pos / REVERSE_TABLE_BUCKET_LENGTH;
					int end_reverse_table_index = (1+section_end_pos) / REVERSE_TABLE_BUCKET_LENGTH;

					start_reverse_table_index = min(start_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH);
					end_reverse_table_index = min(end_reverse_table_index, this_chro_info-> chro_possible_length / REVERSE_TABLE_BUCKET_LENGTH+ 1);

					while(start_reverse_table_index<=end_reverse_table_index)
					{
						search_start = this_chro_info -> reverse_table_start_index [start_reverse_table_index];
						if(search_start<0xffffff00)break;
						start_reverse_table_index++;
					}
					if(search_start>0xffffff00) continue;

					//search_start = this_chro_info -> chro_block_table_start;

					search_end = this_chro_info -> chro_block_table_end;//reverse_table_end_index [end_reverse_table_index];
		
					for(search_block_id=search_start;search_block_id<search_end;search_block_id++){
						if (global_context -> exontable_block_min_start[search_block_id] > section_end_pos) break;
						if (global_context -> exontable_block_max_end[search_block_id] < section_begin_pos) continue;

						int search_item_start = 0, search_item_end = global_context -> exontable_block_end_index[search_block_id];
						if(search_block_id>0)search_item_start = global_context -> exontable_block_end_index[search_block_id-1];

						for(search_item_id = search_item_start ; search_item_id < search_item_end; search_item_id++)
						{
							if (global_context -> exontable_stop[search_item_id] >= section_begin_pos)
							{
								if (global_context -> exontable_start[search_item_id] > section_end_pos) break;
								// there is an overlap >=1 between read and feature.
								// the overlap length is min(end_r, end_F) - max(start_r, start_F) + 1
								
								int is_strand_ok =1;

								if(global_context->is_strand_checked){
									if(global_context->is_strand_checked == 1)
										is_strand_ok = (is_fragment_negative_strand == global_context -> exontable_strand[search_item_id]);
									else// if(global_context->is_strand_checked == 2)
										is_strand_ok = (is_fragment_negative_strand != global_context -> exontable_strand[search_item_id]);
									//SUBREADprintf("%d = %d == %d\n", is_strand_ok, is_fragment_negative_strand, global_context -> exontable_strand[search_item_id]);
								}

								if(is_strand_ok){
									if(nhits<=MAX_HIT_NUMBER - 1)
									{
										hits_indices[nhits] = search_item_id;

										if(global_context -> overlap_length_required !=1)
										{
											int section_overlapped = min(global_context -> exontable_stop[search_item_id] , section_end_pos) 
												- max(global_context -> exontable_start[search_item_id] , section_begin_pos) + 1;
											hits_total_length[nhits] = (short)section_overlapped;
										}

										nhits++;
									}
									else break;
								}
							} 
						}
					}
				}
			}
			else{
				if(global_context->is_split_alignments_only)
				{
					skipped_for_exonic ++;
					if((is_second_read && skipped_for_exonic == 2) || (!global_context -> is_paired_end_data) || (alignment_masks & 0x8))
					{
						if(global_context -> SAM_output_fp)
							fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Nonjunction\t*\t*\n", read_name);

						thread_context->read_counters.unassigned_nonjunction ++;
						return;
					}
				}
			}

			if(is_second_read) nhits2 = nhits;
			else	nhits1 = nhits;
		}
	}	// loop for is_second_read


	if(global_context -> overlap_length_required !=1)
	{
		// merge feature : if a read overlaps with a feature twice or more times (by >=2 segments), the total length of the overlapped bases is calculated and the features with insufficient overlapped bases are removed.
		// both meta-feature mode and feature mode use the same strategy.
		// two ends in a fragment is considered individually; the overlapping bases are not added up.
		int ends;
		for(ends = 0; ends < global_context -> is_paired_end_data + 1; ends++)
		{
			long * hits_indices = ends?hits_indices2:hits_indices1;
			short * hits_total_length = ends?hits_total_length2:hits_total_length1;
			int  nhits = ends?nhits2:nhits1;
			int x1;

			// calculating the summed lengths of overlapping exons
			for(x1=0; x1<nhits;x1++)
			{
				long exon_no = hits_indices[x1];
				if(exon_no>=0x7fffffff) continue;	//already removed
				int x2;
				for(x2=x1+1; x2<nhits; x2++)
				{
					if(hits_indices[x2]==exon_no)
					{
						hits_total_length[x1]+=hits_total_length[x2];
						hits_indices[x2]=0x7fffffff;
					}
				}

				if( global_context -> debug_command [0] == 'D')
					SUBREADprintf("HIT[%d]: OVERLAP=%d\n", x1, hits_total_length[x1]);
				if(hits_total_length[x1]< global_context -> overlap_length_required)
					hits_indices[x1]=0x7fffffff;
			}

			// remove the exons in the hits table when it is marked as removed (0x7fffffff)
			int new_hits=0;
			for(x1 = 0; x1< nhits; x1++)
			{
				if(hits_indices[x1]>=0x7fffffff) continue;
				if(new_hits != x1)
				{
					hits_indices[new_hits]=hits_indices[x1];
					hits_total_length[new_hits]=hits_total_length[x1];		
				}
				new_hits++;
			}
			if(ends) nhits2 = new_hits;
			else nhits1 = new_hits;
		}
	}

	if(nhits2+nhits1==1)
	{
		long hit_exon_id = nhits2?hits_indices2[0]:hits_indices1[0];
		thread_context->count_table[hit_exon_id]++;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;
	}
	else if(nhits2 == 1 && nhits1 == 1 && hits_indices2[0]==hits_indices1[0])
	{
		long hit_exon_id = hits_indices1[0];
		thread_context->count_table[hit_exon_id]++;
		thread_context->nreads_mapped_to_exon++;
		if(global_context -> SAM_output_fp)
		{
			int final_gene_number = global_context -> exontable_geneid[hit_exon_id];
			unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
			fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
		}
		thread_context->read_counters.assigned_reads ++;
	}
	else
	{

		if(nhits2+nhits1>=MAX_HIT_NUMBER)
		{
			print_in_box(80,0,0,"A %s overlapped with %d features.", global_context -> is_paired_end_data?"fragment":"read", nhits2+nhits1);
			print_in_box(80,0,0,"Please increase MAX_HIT_NUMBER in the program");
		}

		long * decision_table_ids = thread_context -> decision_table_ids;
		unsigned char * decision_table_votes = thread_context ->decision_table_votes;
		long * decision_table_exon_ids = thread_context -> decision_table_exon_ids;
		int decision_table_items = 0, xk1, xk2;

		for(is_second_read = 0; is_second_read < 2; is_second_read++)
		{
			if(is_second_read && !global_context -> is_paired_end_data) break;
			long * hits_indices = is_second_read?hits_indices2:hits_indices1;
			int nhits = is_second_read?nhits2:nhits1;
			if (nhits<1) continue;
			if(global_context -> is_gene_level)
			{
				long * uniq_gene_table = thread_context -> uniq_gene_table;
				long * uniq_gene_exonid_table = thread_context -> uniq_gene_exonid_table;
				int uniq_genes = 0;
				for(xk1=0;xk1<nhits;xk1++)
				{
					int gene_id = global_context -> exontable_geneid[hits_indices[xk1]];
					int is_unique = 1;
					for(xk2=0; xk2<uniq_genes; xk2++)
					{
						if(gene_id == uniq_gene_table[xk2])
						{
							is_unique = 0;
							break;
						}
					}
					if(is_unique){
						uniq_gene_exonid_table[uniq_genes] = hits_indices[xk1];
						uniq_gene_table[uniq_genes++] = gene_id;
						if(uniq_genes >= MAX_HIT_NUMBER) break;
					}
				}

				for(xk1=0;xk1<uniq_genes; xk1++)
				{
					long gene_id = uniq_gene_table[xk1];
					int is_fresh = 1;
					if(decision_table_items >= MAX_HIT_NUMBER) break;
					for(xk2=0; xk2<decision_table_items; xk2++)
					{
						if(gene_id == decision_table_ids[xk2])
						{
							decision_table_votes[xk2]++;
							is_fresh = 0;
							break;
						}
						
					}
					if(is_fresh)
					{
						decision_table_votes[decision_table_items] = 1;
						decision_table_exon_ids[decision_table_items] = uniq_gene_exonid_table[xk1];
						decision_table_ids[decision_table_items++] = gene_id;
					}
				}
			}
			else
			{
				for(xk1=0;xk1<nhits;xk1++)
				{
					long exon_id = hits_indices[xk1];
					//long exon_id = global_context -> exontable_geneid[hits_indices[xk1]];
					int is_fresh = 1;
					if(decision_table_items >= MAX_HIT_NUMBER) break;
					for(xk2=0; xk2<decision_table_items; xk2++)
					{
						if(exon_id == decision_table_ids[xk2])
						{
							decision_table_votes[xk2]++;
							is_fresh = 0;
							break;
						}
					}
					if(is_fresh)
					{
						decision_table_votes[decision_table_items] = 1;
						decision_table_ids[decision_table_items++] = exon_id;
					}

				}
			}

		}

		/*
		printf("R has %d hits\n", decision_table_items);
		int i;
		for(i=0; i<decision_table_items; i++)
		{
			long geneid = global_context -> exontable_geneid[decision_table_ids[i]] ;
			printf("%ld : %ld,   ",decision_table_ids[i], geneid);
		}
		puts("");

		*/

		if(decision_table_items>0)
		{
			int max_votes = 0;
			int top_voters = 0;
			long top_voter_id = 0;

			for(xk1 = 0; xk1 < decision_table_items; xk1++)
			{
				if(decision_table_votes[xk1] > max_votes)
				{
					max_votes = decision_table_votes[xk1];
					top_voters = 1;
					top_voter_id = (global_context -> is_gene_level)?decision_table_exon_ids[xk1]:decision_table_ids[xk1];
				}
				else
					if(decision_table_votes[xk1] == max_votes) top_voters++;
			}

			if(top_voters == 1 && (!global_context -> is_multi_overlap_allowed))
			{
				thread_context->count_table[top_voter_id]++;
				thread_context->nreads_mapped_to_exon++;
				if(global_context -> SAM_output_fp)
				{
					int final_gene_number = global_context -> exontable_geneid[top_voter_id];
					unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
					if(decision_table_items>1)
						// assigned by read-voting
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tVotes/Targets=%d/%d\n", read_name, final_feture_name, max_votes, decision_table_items);
					else
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\t*\n", read_name, final_feture_name);
				}
				thread_context->read_counters.assigned_reads ++;

			}
			else if(top_voters >1 || (global_context -> is_multi_overlap_allowed))
			{
				if(global_context -> is_multi_overlap_allowed)
				{
					char final_feture_names[1000];
					int assigned_no = 0;
					final_feture_names[0]=0;
					for(xk1 = 0; xk1 < decision_table_items; xk1++)
					{
						//if(decision_table_votes[xk1] == max_votes)
						if(decision_table_votes[xk1] >= 1)
						{
							long tmp_voter_id = (global_context -> is_gene_level)?decision_table_exon_ids[xk1]:decision_table_ids[xk1];
							//printf("TVID=%ld; HITID=%d\n", tmp_voter_id, xk1);
							thread_context->count_table[tmp_voter_id]++;

							if(global_context -> SAM_output_fp)
							{
								if(strlen(final_feture_names)<700) 
								{
										int final_gene_number = global_context -> exontable_geneid[tmp_voter_id];
										unsigned char * final_feture_name = global_context -> gene_name_array[final_gene_number];
										strncat(final_feture_names, (char *)final_feture_name, 999);
										strncat(final_feture_names, ",", 999);
										assigned_no++;
								}
							}
						}
					}
					final_feture_names[999]=0;
					thread_context->nreads_mapped_to_exon++;
					if(global_context -> SAM_output_fp)
					{
						int ffnn = strlen(final_feture_names);
						if(ffnn>0) final_feture_names[ffnn-1]=0;
						// overlapped but still assigned 
						fprintf(global_context -> SAM_output_fp,"%s\tAssigned\t%s\tTotal=%d\n", read_name, final_feture_names, assigned_no);
					}
					thread_context->read_counters.assigned_reads ++;

				}
				else{
					if(global_context -> SAM_output_fp)
						fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_Ambiguit\t*\tNumber_Of_Overlapped_Genes=%d\n", read_name, top_voters);

					thread_context->read_counters.unassigned_ambiguous ++;
				}
			}
		}
		else{
			if(global_context -> SAM_output_fp)
				fprintf(global_context -> SAM_output_fp,"%s\tUnassigned_NoFeatures\t*\t*\n", read_name);

			thread_context->read_counters.unassigned_nofeatures ++;
		}
	}

}

void * feature_count_worker(void * vargs)
{
	void ** args = (void **) vargs;

	fc_thread_global_context_t * global_context = args[0];
	fc_thread_thread_context_t * thread_context = args[1];

	free(vargs);


	//printf("QQQ0:T%d\n", thread_context->thread_id);
	//Rprintf("QQQ1:T%d\n", thread_context->thread_id);
	//printf("QQQ2:T%d\n", thread_context->thread_id);

	if(global_context -> is_SAM_file)
	{
		//thread_context -> current_read_length1 = global_context -> read_length;
		//thread_context -> current_read_length2 = global_context -> read_length;
		while (1)
		{
			while(1)
			{
				int is_retrieved = 0;
				pthread_spin_lock(&thread_context->input_buffer_lock);
				if(thread_context->input_buffer_remainder)
				{
					int is_second_read;
					unsigned int buffer_read_bytes ;
					unsigned int buffer_read_ptr;
					if(thread_context->input_buffer_remainder <= thread_context->input_buffer_write_ptr)
						buffer_read_ptr = thread_context->input_buffer_write_ptr - thread_context->input_buffer_remainder; 
					else
						buffer_read_ptr = thread_context->input_buffer_write_ptr + global_context->input_buffer_max_size - thread_context->input_buffer_remainder;

					//if(buffer_read_ptr>= global_context->input_buffer_max_size)
					//	if(buffer_read_ptr>6*1024*1024) printf("REALLY BIG PTR:%u = %u + %u - %u\n", buffer_read_ptr, thread_context->input_buffer_write_ptr , global_context->input_buffer_max_size, thread_context->input_buffer_remainder);

					for(is_second_read = 0; is_second_read < (global_context->is_paired_end_data ? 2:1); is_second_read++)
					{
						char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;
						//printf("R=%u; WPTR=%u ;RPTR=%u\n", thread_context->input_buffer_remainder, thread_context->input_buffer_write_ptr, buffer_read_ptr);
						//if(buffer_read_ptr % 7 == 0)
						//	fflush(stdout);
						
						for(buffer_read_bytes=0; ; buffer_read_bytes++)
						{
							//printf("%p + %d\n", thread_context->input_buffer, buffer_read_ptr);
							//if(buffer_read_ptr>6*1024*1024) printf("VERY BIG PTR:%u > %u\n", buffer_read_ptr , global_context->input_buffer_max_size);
							char nch =  thread_context->input_buffer[buffer_read_ptr ++];
							curr_line_buff[buffer_read_bytes] = nch;
							if(buffer_read_ptr >= global_context->input_buffer_max_size)
								buffer_read_ptr = 0; 
							if(nch=='\n' || buffer_read_bytes>2998){
								curr_line_buff[buffer_read_bytes+1]=0;
								break;
							}
						}

						//printf("%s\n", curr_line_buff);

						//if(buffer_read_bytes + 1 > thread_context->input_buffer_remainder)
						//	(*(int*)NULL) = 1;
						thread_context->input_buffer_remainder -= buffer_read_bytes + 1;
					}
					is_retrieved = 1;

				}

				pthread_spin_unlock(&thread_context->input_buffer_lock);
				if(global_context->is_all_finished && !is_retrieved) return NULL;

				if(is_retrieved) break;
				else
					usleep(tick_time);
			}



			process_line_buffer(global_context, thread_context);

		}
	}
	else
	{	// if is BAM: decompress the chunk and process reads.
		char * PDATA = malloc(2*70000);
		SamBam_Alignment * aln = &thread_context->aln_buffer;

		//thread_context -> current_read_length1 = global_context -> read_length;
		//thread_context -> current_read_length2 = global_context -> read_length;
		while(1)
		{
			int PDATA_len = 0;
			while(1)
			{
				int is_retrieved = 0;
				PDATA_len = 0;
				//retrieve the next chunk.

				pthread_spin_lock(&thread_context->input_buffer_lock);
				if(thread_context->input_buffer_remainder)
				{
					assert(thread_context->input_buffer_remainder>4);
					unsigned int tail_bytes = global_context->input_buffer_max_size - thread_context -> chunk_read_ptr ;
					if(tail_bytes<4)
					{
						thread_context -> chunk_read_ptr = 0;
						thread_context -> input_buffer_remainder -= tail_bytes;
						memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
					}
					else
					{
						memcpy(&PDATA_len, thread_context->input_buffer + thread_context -> chunk_read_ptr , 4);
						if(PDATA_len==0)
						{
							thread_context -> chunk_read_ptr = 0;
							thread_context -> input_buffer_remainder -= tail_bytes;
							memcpy(&PDATA_len, thread_context->input_buffer , 4);
						}
					}
					thread_context -> chunk_read_ptr+=4;
					thread_context -> input_buffer_remainder -= 4;

					//fprintf(stderr,"chunk_read_ptr=%d , input_buffer_remainder = %d\n", thread_context -> chunk_read_ptr , thread_context -> input_buffer_remainder);
					if(PDATA_len<0 || PDATA_len > 140000)
					{
						SUBREADprintf("THREAD ABNORMALLY QUIT\n");
						return NULL;
					}

					memcpy(PDATA, thread_context -> input_buffer + thread_context -> chunk_read_ptr , PDATA_len);
					thread_context -> chunk_read_ptr += PDATA_len;
					thread_context -> input_buffer_remainder -= PDATA_len;

					if( PDATA_len > 0 )
						is_retrieved = 1;
				}


				pthread_spin_unlock(&thread_context->input_buffer_lock);
				if(global_context->is_all_finished && !is_retrieved){
					free(PDATA);
					return NULL;
				}

				if(is_retrieved) break;
				else
					usleep(tick_time);

			}
			
			// convert binary reads into sam lines and process;
			int processed_reads = 0, PDATA_ptr = 0;
			while(PDATA_ptr < PDATA_len)
			{
				int is_second_read;
				for(is_second_read = 0; is_second_read <= global_context -> is_paired_end_data; is_second_read++)
				{
					int binary_read_len, local_PDATA_ptr = PDATA_ptr;
					char * curr_line_buff = is_second_read?thread_context -> line_buffer2:thread_context -> line_buffer1;

					memcpy(&binary_read_len, PDATA + PDATA_ptr, 4);
					int ret = PBam_chunk_gets(PDATA, &local_PDATA_ptr, PDATA_len, global_context -> sambam_chro_table, curr_line_buff, 2999, aln,0);
					//printf("LL=%s\n", curr_line_buff);
					if(ret<0)
						SUBREADprintf("READ DECODING ERROR!\n");

					PDATA_ptr += 4+binary_read_len;
					processed_reads++;
				}

				process_line_buffer(global_context, thread_context);
				//printf("LE\n\n");
			}
		}
	}
}

void fc_thread_merge_results(fc_thread_global_context_t * global_context, unsigned int * nreads , unsigned long long int *nreads_mapped_to_exon, fc_read_counters * my_read_counter)
{
	int xk1, xk2;

	long long int total_input_reads = 0 ;

	(*nreads_mapped_to_exon)=0;

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		for(xk2=0; xk2<global_context -> exontable_exons; xk2++)
		{
			nreads[xk2]+=global_context -> thread_contexts[xk1].count_table[xk2];
		}
		total_input_reads += global_context -> thread_contexts[xk1].all_reads;
		(*nreads_mapped_to_exon) += global_context -> thread_contexts[xk1].nreads_mapped_to_exon;

		global_context -> read_counters.unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		global_context -> read_counters.unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		global_context -> read_counters.unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		global_context -> read_counters.unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		global_context -> read_counters.unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		global_context -> read_counters.unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		global_context -> read_counters.unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		global_context -> read_counters.unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		global_context -> read_counters.unassigned_nonjunction += global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction;
		global_context -> read_counters.unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		global_context -> read_counters.assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;

		my_read_counter->unassigned_ambiguous += global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous;
		my_read_counter->unassigned_nofeatures += global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures;
		my_read_counter->unassigned_unmapped += global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped;
		my_read_counter->unassigned_mappingquality += global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality;
		my_read_counter->unassigned_fragmentlength += global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength;
		my_read_counter->unassigned_chimericreads += global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads;
		my_read_counter->unassigned_multimapping += global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping;
		my_read_counter->unassigned_secondary += global_context -> thread_contexts[xk1].read_counters.unassigned_secondary;
		my_read_counter->unassigned_nonjunction += global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction;
		my_read_counter->unassigned_duplicate += global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate;
		my_read_counter->assigned_reads += global_context -> thread_contexts[xk1].read_counters.assigned_reads;
		
	}

	char pct_str[10];
	if(total_input_reads>0)
		sprintf(pct_str,"(%.1f%%%%)", (*nreads_mapped_to_exon)*100./total_input_reads);
	else	pct_str[0]=0;

	print_in_box(80,0,0,"   Total %s : %llu", global_context -> is_paired_end_data?"fragments":"reads", total_input_reads); 
	print_in_box(pct_str[0]?81:80,0,0,"   Successfully assigned %s : %llu %s", global_context -> is_paired_end_data?"fragments":"reads", *nreads_mapped_to_exon,pct_str); 
	print_in_box(80,0,0,"   Running time : %.2f minutes", (miltime() - global_context -> start_time)/60);
	print_in_box(80,0,0,"");
}

HashTable * load_alias_table(char * fname)
{
	FILE * fp = f_subr_open(fname, "r");
	if(!fp)
	{
		print_in_box(80,0,0,"WARNING unable to open alias file '%s'", fname);
		return NULL;
	}

	char * fl = malloc(2000);

	HashTable * ret = HashTableCreate(1013);
	HashTableSetDeallocationFunctions(ret, free, free);
	HashTableSetKeyComparisonFunction(ret, fc_strcmp);
	HashTableSetHashFunction(ret, fc_chro_hash);
	
	while (1)
	{
		char *ret_fl = fgets(fl, 1999, fp);
		if(!ret_fl) break;
		if(fl[0]=='#') continue;
		char * sam_chr = NULL;
		char * anno_chr = strtok_r(fl, ",", &sam_chr);
		if((!sam_chr)||(!anno_chr)) continue;

		sam_chr[strlen(sam_chr)-1]=0;
		char * anno_chr_buf = malloc(strlen(anno_chr)+1);
		strcpy(anno_chr_buf, anno_chr);
		char * sam_chr_buf = malloc(strlen(sam_chr)+1);
		strcpy(sam_chr_buf, sam_chr);
		
		//printf("ALIAS: %s -> %s\n", sam_chr, anno_chr);
		HashTablePut(ret, sam_chr_buf, anno_chr_buf);
	}


	free(fl);
	return ret;
}

void fc_thread_init_global_context(fc_thread_global_context_t * global_context, unsigned int buffer_size, unsigned short threads, int line_length , int is_PE_data, int min_pe_dist, int max_pe_dist, int is_gene_level, int is_overlap_allowed, int is_strand_checked, char * output_fname, int is_sam_out, int is_both_end_required, int is_chimertc_disallowed, int is_PE_distance_checked, char *feature_name_column, char * gene_id_column, int min_map_qual_score, int is_multi_mapping_allowed, int is_SAM, char * alias_file_name, char * cmd_rebuilt, int is_input_file_resort_needed, int feature_block_size, int isCVersion, int fiveEndExtension,  int threeEndExtension, int minReadOverlap, int is_split_alignments_only, int reduce_5_3_ends_to_one, char * debug_command, int is_duplicate_ignored)
{

	global_context -> input_buffer_max_size = buffer_size;
	global_context -> all_reads = 0;
	global_context -> redo = 0;
	global_context -> SAM_output_fp = NULL;


	global_context -> isCVersion = isCVersion;
	global_context -> is_read_details_out = is_sam_out;
	global_context -> is_multi_overlap_allowed = is_overlap_allowed;
	global_context -> is_paired_end_data = is_PE_data;
	global_context -> is_gene_level = is_gene_level;
	global_context -> is_strand_checked = is_strand_checked;
	global_context -> is_both_end_required = is_both_end_required;
	global_context -> is_chimertc_disallowed = is_chimertc_disallowed;
	global_context -> is_PE_distance_checked = is_PE_distance_checked;
	global_context -> is_multi_mapping_allowed = is_multi_mapping_allowed;
	global_context -> is_split_alignments_only = is_split_alignments_only;
	global_context -> is_duplicate_ignored = is_duplicate_ignored;
	global_context -> reduce_5_3_ends_to_one = reduce_5_3_ends_to_one;
	global_context -> is_SAM_file = is_SAM;


	global_context -> thread_number = threads;
	global_context -> min_mapping_quality_score = min_map_qual_score;
	global_context -> unistr_buffer_size = 1024*1024*2;
	global_context -> unistr_buffer_used = 0;
	global_context -> unistr_buffer_space = malloc(global_context -> unistr_buffer_size);
	global_context -> annot_chro_name_alias_table = NULL;
	global_context -> cmd_rebuilt = cmd_rebuilt;
	global_context -> is_input_file_resort_needed = is_input_file_resort_needed;
	global_context -> feature_block_size = feature_block_size;
	global_context -> five_end_extension = fiveEndExtension;
	global_context -> three_end_extension = threeEndExtension;
	global_context -> overlap_length_required = minReadOverlap;
	global_context -> debug_command = debug_command;

	global_context -> read_counters.unassigned_ambiguous=0;
	global_context -> read_counters.unassigned_nofeatures=0;
	global_context -> read_counters.unassigned_unmapped=0;
	global_context -> read_counters.unassigned_mappingquality=0;
	global_context -> read_counters.unassigned_fragmentlength=0;
	global_context -> read_counters.unassigned_chimericreads=0;
	global_context -> read_counters.unassigned_multimapping=0;
	global_context -> read_counters.unassigned_secondary=0;
	global_context -> read_counters.unassigned_nonjunction=0;
	global_context -> read_counters.unassigned_duplicate=0;
	global_context -> read_counters.assigned_reads=0;

	if(alias_file_name && alias_file_name[0])
	{
		strcpy(global_context -> alias_file_name,alias_file_name);
		global_context -> annot_chro_name_alias_table = load_alias_table(alias_file_name);
	}
	else	global_context -> alias_file_name[0]=0;

	strcpy(global_context -> feature_name_column,feature_name_column);
	strcpy(global_context -> gene_id_column,gene_id_column);
	strcpy(global_context -> output_file_name, output_fname);

	global_context -> min_paired_end_distance = min_pe_dist;
	global_context -> max_paired_end_distance = max_pe_dist;
	global_context -> thread_number = threads;
	global_context -> line_length = line_length;


}
int fc_thread_start_threads(fc_thread_global_context_t * global_context, int et_exons, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, char * et_anno_chr_2ch, char ** et_anno_chrs, long * et_anno_chr_heads, long * et_bk_end_index, long * et_bk_min_start, long * et_bk_max_end, int read_length)
{
	int xk1;

	global_context -> read_length = read_length;
	global_context -> is_unpaired_warning_shown = 0;
	global_context -> is_stake_warning_shown = 0;


	if(global_context -> is_read_details_out)
	{
		char tmp_fname[350];
		sprintf(tmp_fname, "%s.featureCounts", global_context -> raw_input_file_name);
		global_context -> SAM_output_fp = f_subr_open(tmp_fname, "w");
		if(!global_context -> SAM_output_fp)
		{
			SUBREADprintf("Unable to create file '%s'; the read assignment details are not written.\n", tmp_fname);
		}
	}
	else
		global_context -> SAM_output_fp = NULL;

	global_context -> redo = 0;
	global_context -> exontable_geneid = et_geneid;
	global_context -> exontable_chr = et_chr;
	global_context -> exontable_start = et_start;
	global_context -> exontable_stop = et_stop;
	global_context -> exontable_strand = (char *)et_strand;
	global_context -> exontable_anno_chr_2ch = et_anno_chr_2ch;
	global_context -> exontable_anno_chrs = et_anno_chrs;
	global_context -> exontable_anno_chr_heads = et_anno_chr_heads;
	global_context -> exontable_block_end_index = et_bk_end_index;
	global_context -> exontable_block_max_end = et_bk_max_end;
	global_context -> exontable_block_min_start = et_bk_min_start;

	global_context -> is_all_finished = 0;
	global_context -> thread_contexts = malloc(sizeof(fc_thread_thread_context_t) * global_context -> thread_number);
	for(xk1=0; xk1<global_context -> thread_number; xk1++)
	{
	//	printf("CHRR_MALLOC\n");
		pthread_spin_init(&global_context->thread_contexts[xk1].input_buffer_lock, PTHREAD_PROCESS_PRIVATE);
		global_context -> thread_contexts[xk1].input_buffer_remainder = 0;
		global_context -> thread_contexts[xk1].input_buffer_write_ptr = 0;
		global_context -> thread_contexts[xk1].input_buffer = malloc(global_context -> input_buffer_max_size);
		global_context -> thread_contexts[xk1].thread_id = xk1;
		global_context -> thread_contexts[xk1].chunk_read_ptr = 0;
		global_context -> thread_contexts[xk1].count_table = calloc(sizeof(unsigned int), et_exons);
		global_context -> thread_contexts[xk1].nreads_mapped_to_exon = 0;
		global_context -> thread_contexts[xk1].all_reads = 0;
		global_context -> thread_contexts[xk1].line_buffer1 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].line_buffer2 = malloc(global_context -> line_length + 2);
		global_context -> thread_contexts[xk1].chro_name_buff = malloc(CHROMOSOME_NAME_LENGTH);
		global_context -> thread_contexts[xk1].strm_buffer = malloc(sizeof(z_stream));

		global_context -> thread_contexts[xk1].read_counters.assigned_reads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_ambiguous = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nofeatures = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_unmapped = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_mappingquality = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_fragmentlength = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_chimericreads = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_multimapping = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_secondary = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_nonjunction = 0;
		global_context -> thread_contexts[xk1].read_counters.unassigned_duplicate = 0;

		if(!global_context ->  thread_contexts[xk1].count_table) return 1;
		void ** thread_args = malloc(sizeof(void *)*2);
		thread_args[0] = global_context;
		thread_args[1] = & global_context -> thread_contexts[xk1];

		if(global_context -> thread_number>1 || ! global_context -> is_SAM_file)
			pthread_create(&global_context -> thread_contexts[xk1].thread_object, NULL, feature_count_worker, thread_args);
	}

	return 0;
}

void fc_thread_destroy_thread_context(fc_thread_global_context_t * global_context)
{
	int xk1;
	if(global_context -> is_read_details_out)
	{
		fclose(global_context -> SAM_output_fp);
		global_context -> SAM_output_fp = NULL;
	}

	for(xk1=0; xk1<global_context-> thread_number; xk1++)
	{
		//printf("CHRR_FREE\n");
		free(global_context -> thread_contexts[xk1].count_table);	
		free(global_context -> thread_contexts[xk1].line_buffer1);	
		free(global_context -> thread_contexts[xk1].line_buffer2);	
		free(global_context -> thread_contexts[xk1].input_buffer);
		free(global_context -> thread_contexts[xk1].chro_name_buff);
		free(global_context -> thread_contexts[xk1].strm_buffer);
		pthread_spin_destroy(&global_context -> thread_contexts[xk1].input_buffer_lock);
	}
	free(global_context -> thread_contexts);
}
void fc_thread_wait_threads(fc_thread_global_context_t * global_context)
{
	int xk1;
	for(xk1=0; xk1<global_context-> thread_number; xk1++)
		pthread_join(global_context -> thread_contexts[xk1].thread_object, NULL);
}

int resort_input_file(fc_thread_global_context_t * global_context)
{
	char * temp_file_name = malloc(300), * fline = malloc(3000);
	SamBam_FILE * sambam_reader ;

	if(!global_context->redo)
		print_in_box(80,0,0,"   Resort the input file ...");
	sprintf(temp_file_name, "./temp-core-%06u-%08X.sam", getpid(), rand());
	sambam_reader = SamBam_fopen(global_context-> input_file_name, global_context-> is_SAM_file?SAMBAM_FILE_SAM:SAMBAM_FILE_BAM);

	if(!sambam_reader){
		SUBREADprintf("Unable to open %s.\n", global_context-> input_file_name);
		return -1;
	}

	SAM_sort_writer writer;
	int ret = sort_SAM_create(&writer, temp_file_name, ".");
	if(ret)
	{
		SUBREADprintf("Unable to sort input file because temporary file '%s' cannot be created.\n", temp_file_name);
		return -1;
	}
	int is_read_len_warned = 0;

	while(1)
	{
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) break;
		int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
		if(ret<0) 
		{
			if(!is_read_len_warned)
				print_in_box(80,0,0,"WARNING: reads with very long names were found.");
			is_read_len_warned = 1;
		//	break;
		}
	//printf("N1=%llu\n",  writer.unpaired_reads);
	}

	sort_SAM_finalise(&writer);
	print_in_box(80,0,0,"   %llu read%s ha%s missing mates.", writer.unpaired_reads, writer.unpaired_reads>1?"s":"", writer.unpaired_reads>1?"ve":"s");
	print_in_box(80,0,0,"   Input was converted to a format accepted by featureCounts.");

	SamBam_fclose(sambam_reader);
	strcpy(global_context-> input_file_name, temp_file_name);
	global_context->is_SAM_file = 1;
	free(temp_file_name);
	free(fline);
	return 0;
}


void fc_write_final_gene_results(fc_thread_global_context_t * global_context, int * et_geneid, char ** et_chr, long * et_start, long * et_stop, unsigned char * et_strand, const char * out_file, int features, unsigned int ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	int xk1;
	int genes = global_context -> gene_name_table -> numOfElements;
	unsigned int *gene_columns;

	FILE * fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
		return;
	}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}

	char * tmp_ptr = NULL, * next_fn;
	int non_empty_files = 0, i_files=0;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
		{
			fprintf(fp_out,"\t%s", next_fn);
			non_empty_files ++;
		}
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");

	gene_columns = calloc(sizeof(unsigned int) , genes * non_empty_files);
	unsigned int * gene_exons_number = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_pointer = calloc(sizeof(unsigned int) , genes);
	unsigned int * gene_exons_start = malloc(sizeof(unsigned int) * features);
	unsigned int * gene_exons_end = malloc(sizeof(unsigned int) * features);
	char ** gene_exons_chr = malloc(sizeof(char *) * features);
	char * gene_exons_strand = malloc(features);

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		gene_exons_number[gene_id]++;
	}

	unsigned int accumulative_no = 0;
	unsigned longest_gene_exons = 0;
	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		unsigned int tmpv = gene_exons_number[xk1];
		longest_gene_exons = max(longest_gene_exons, tmpv);
		accumulative_no += gene_exons_number[xk1];
		gene_exons_number[xk1] = accumulative_no - tmpv;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1];
		int gene_write_ptr = gene_exons_number[gene_id] + gene_exons_pointer[gene_id];

		gene_exons_chr[gene_write_ptr] = et_chr[xk1];
		gene_exons_start[gene_write_ptr] = et_start[xk1]; 
		gene_exons_end[gene_write_ptr] = et_stop[xk1]; 
		gene_exons_strand[gene_write_ptr] = et_strand[xk1]; 

		gene_exons_pointer[gene_id]++;
	}

	for(xk1 = 0; xk1 < features; xk1++)
	{
		int gene_id = et_geneid[xk1], k_noempty = 0;
		for(i_files=0;i_files < n_input_files; i_files++)
		{
			if(column_numbers[i_files]==NULL) continue;
			gene_columns[gene_id * non_empty_files + k_noempty ] += column_numbers[i_files][xk1];
			k_noempty++;
		}
	}


	char *is_occupied = malloc(longest_gene_exons);
	unsigned int * input_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);
	unsigned int * output_start_stop_list = malloc(longest_gene_exons * sizeof(int) * 2);

	char * out_chr_list = malloc(longest_gene_exons * (1+global_context -> longest_chro_name) + 1);
	char * out_start_list = malloc(11 * longest_gene_exons + 1);
	char * out_end_list = malloc(11 * longest_gene_exons + 1);
	char * out_strand_list = malloc(2 * longest_gene_exons + 1);

	for(xk1 = 0 ; xk1 < genes; xk1++)
	{
		int xk2;
		
		memset(is_occupied,0,gene_exons_pointer[xk1]);
		out_chr_list[0]=0;
		out_start_list[0]=0;
		out_end_list[0]=0;
		out_strand_list[0]=0;
		int gene_nonoverlap_len =0;

		for(xk2=0; xk2<gene_exons_pointer[xk1]; xk2++)
		{
			if(!is_occupied[xk2])
			{
				int xk3;
				char * matched_chr = gene_exons_chr[xk2 + gene_exons_number[xk1]];
				char matched_strand = gene_exons_strand[xk2 + gene_exons_number[xk1]];

				memset(input_start_stop_list, 0, gene_exons_pointer[xk1] * sizeof(int) * 2);
				int gap_merge_ptr = 1;
				input_start_stop_list[0] = gene_exons_start[xk2 + gene_exons_number[xk1]];
				input_start_stop_list[1] = gene_exons_end[xk2 + gene_exons_number[xk1]] + 1;

				for(xk3 = xk2+1; xk3 < gene_exons_pointer[xk1]; xk3++)
				{
					if((!is_occupied[xk3]) && strcmp(matched_chr, gene_exons_chr[xk3+gene_exons_number[xk1]])==0 && matched_strand == gene_exons_strand[xk3 + gene_exons_number[xk1]])
					{
						is_occupied[xk3]=1;
						input_start_stop_list[gap_merge_ptr*2] = gene_exons_start[xk3+gene_exons_number[xk1]]; 
						input_start_stop_list[gap_merge_ptr*2+1] = gene_exons_end[xk3+gene_exons_number[xk1]]+1;
						gap_merge_ptr++;
					}
				}

				{
						int merged_gaps = mergeIntervals(input_start_stop_list, output_start_stop_list, gap_merge_ptr);

						for(xk3=0; xk3<merged_gaps; xk3++)
						{
							char numbbuf[12];
							strcat(out_chr_list, matched_chr);
							strcat(out_chr_list, ";");

							sprintf(numbbuf,"%u;", output_start_stop_list[xk3 * 2]);
							strcat(out_start_list, numbbuf);
							sprintf(numbbuf,"%u;", output_start_stop_list[xk3 * 2 + 1] - 1);
							strcat(out_end_list, numbbuf);
							sprintf(numbbuf,"%c;", matched_strand?'-':'+');
							strcat(out_strand_list, numbbuf);

							gene_nonoverlap_len += output_start_stop_list[xk3 * 2 + 1] - output_start_stop_list[xk3 * 2];
						}
				}	
			}
		}

		unsigned char * gene_symbol = global_context -> gene_name_array [xk1];

		#define _cut_tail(x) (x)[strlen(x)-1]=0

		_cut_tail(out_chr_list);
		_cut_tail(out_start_list);
		_cut_tail(out_end_list);
		_cut_tail(out_strand_list);

		fprintf(fp_out, "%s\t%s\t%s\t%s\t%s\t%d"    , gene_symbol, out_chr_list, out_start_list, out_end_list, out_strand_list, gene_nonoverlap_len);

		// all exons: gene_exons_number[xk1] : gene_exons_pointer[xk1]
		for(i_files=0; i_files< non_empty_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				fprintf(fp_out,"\t%u", gene_columns[i_files+non_empty_files*xk1]);
			}
		}
		fprintf(fp_out,"\n");

	}
	free(is_occupied);
	free(input_start_stop_list);
	free(output_start_stop_list);
	free(out_chr_list);
	free(out_strand_list);
	free(out_start_list);
	free(out_end_list);

	free(gene_exons_number);
	free(gene_exons_pointer);
	free(gene_columns);
	free(gene_exons_chr);
	free(gene_exons_start);
	free(gene_exons_end);
	free(gene_exons_strand);
	fclose(fp_out);
}

void fc_write_final_counts(fc_thread_global_context_t * global_context, const char * out_file, int nfiles, char * file_list, unsigned int ** column_numbers, fc_read_counters *read_counters, int isCVersion)
{
	char fname[300];
	int i_files, xk1;

	sprintf(fname, "%s.summary", out_file);
	FILE * fp_out = f_subr_open(fname,"w");

	if(!fp_out){
		SUBREADprintf("Unable to create summary file '%s'\n", fname);
		return;
	}

	fprintf(fp_out,"Status");
	char * next_fn = file_list;
	
	for(i_files=0; i_files<nfiles; i_files++)
	{
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);

		next_fn += strlen(next_fn)+1;
	}

	fprintf(fp_out,"\n");
	char * keys [] ={ "Assigned" , "Unassigned_Ambiguity", "Unassigned_MultiMapping" ,"Unassigned_NoFeatures", "Unassigned_Unmapped", "Unassigned_MappingQuality", "Unassigned_FragementLength", "Unassigned_Chimera", "Unassigned_Secondary", "Unassigned_Nonjunction", "Unassigned_Duplicate"};

	for(xk1=0; xk1<11; xk1++)
	{
		fprintf(fp_out,"%s", keys[xk1]);
		for(i_files = 0; i_files < nfiles; i_files ++)
		{
			unsigned long long * array_0 = (unsigned long long *)&(read_counters[i_files]);
			unsigned long long * cntr = array_0 + xk1;
			if(column_numbers[i_files])
			{
				fprintf(fp_out,"\t%llu", *(cntr));
			}
		}
		fprintf(fp_out,"\n");
	}


	fclose(fp_out);
}
void fc_write_final_results(fc_thread_global_context_t * global_context, const char * out_file, int features, unsigned int ** column_numbers, char * file_list, int n_input_files, fc_feature_info_t * loaded_features, int header_out)
{
	/* save the results */
	FILE * fp_out;
	int i, i_files = 0;
	fp_out = f_subr_open(out_file,"w");
	if(!fp_out){
		SUBREADprintf("Failed to create file %s\n", out_file);
			return;
		}

	if(header_out)
	{
		fprintf(fp_out, "# Program:featureCounts v%s", SUBREAD_VERSION);
		if(global_context->cmd_rebuilt)
			fprintf(fp_out, "; Command:%s", global_context->cmd_rebuilt);
		fprintf(fp_out, "\n");
	}



	char * tmp_ptr = NULL, * next_fn;
	fprintf(fp_out,"Geneid\tChr\tStart\tEnd\tStrand\tLength");
	next_fn = strtok_r(file_list, ";", &tmp_ptr);
	while(1){
		if(!next_fn||strlen(next_fn)<1) break;
		if(column_numbers[i_files])
			fprintf(fp_out,"\t%s", next_fn);
		next_fn = strtok_r(NULL, ";", &tmp_ptr);
		i_files++;
	}
	fprintf(fp_out,"\n");
	for(i=0;i<features;i++)
	{
		fprintf(fp_out,"%s\t%s\t%u\t%u\t%c\t%d", global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos,
 							   global_context -> unistr_buffer_space + loaded_features[i].feature_name_pos + loaded_features[i].chro_name_pos_delta,
						   	   loaded_features[i].start, loaded_features[i].end, loaded_features[i].is_negative_strand?'-':'+',loaded_features[i].end-loaded_features[i].start+1);
		for(i_files=0; i_files<n_input_files; i_files++)
		{
			if(column_numbers[i_files])
			{
				int sorted_exon_no = loaded_features[i].sorted_order;
				fprintf(fp_out,"\t%d", column_numbers[i_files][sorted_exon_no]);
			}
		}
		fprintf(fp_out,"\n");
	}

	fclose(fp_out);
}

static struct option long_options[] =
{
	{"primary",no_argument, 0, 0},
	{"readExtension5", required_argument, 0, 0},
	{"readExtension3", required_argument, 0, 0},
	{"read2pos", required_argument, 0, 0},
	{"minReadOverlap", required_argument, 0, 0},
	{"countSplitAlignmentsOnly", no_argument, 0, 0},
	{"debugCommand", required_argument, 0, 0},
	{"ignoreDup", no_argument, 0, 0},
	{0, 0, 0, 0}
};

void print_usage()
{
	SUBREADprintf("\nVersion %s\n\n", SUBREAD_VERSION);

	SUBREADputs("\nUsage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... \n");
	SUBREADputs("    Required parameters:\n"); 
	SUBREADputs("    -a <input>\tGive the name of the annotation file. The program assumes"); 
	SUBREADputs("              \tthat the provided annotation file is in GTF format. Use -F"); 
	SUBREADputs("              \toption to specify other annotation formats."); 
	SUBREADputs("    "); 
	SUBREADputs("    -o <input>\tGive the name of the output file. The output file contains"); 
	SUBREADputs("              \tthe number of reads assigned to each meta-feature (or each"); 
	SUBREADputs("              \tfeature if -f is specified). A meta-feature is the aggregation");
	SUBREADputs("              \tof features, grouped by using gene identifiers. Please refer");
	SUBREADputs("              \tto the users guide for more details."); 
	SUBREADputs("    "); 
	SUBREADputs("   input_files\tGive the names of input read files that include the read");
	SUBREADputs("              \tmapping results. Format of input files is automatically");
	SUBREADputs("              \tdetermined (SAM or BAM). Paired-end reads will be");
	SUBREADputs("              \tautomatically re-ordered if it is found that reads from the");
	SUBREADputs("              \tsame pair are not adjacent to each other. Multiple files can");
	SUBREADputs("              \tbe provided at the same time."); 
	SUBREADputs("    "); 
	SUBREADputs("    Optional parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -A <input>\tSpecify the name of a file including aliases of chromosome");
	SUBREADputs("              \tnames. The file should be a comma delimited text file that");
	SUBREADputs("              \tincludes two columns. The first column gives the chromosome");
	SUBREADputs("              \tnames used in the annotation and the second column gives the");
	SUBREADputs("              \tchromosome names used by reads. This file should not contain");
	SUBREADputs("              \theader lines. Names included in this file are case sensitive.");
	SUBREADputs("    "); 
	SUBREADputs("    -F <input>\tSpecify the format of the annotation file. Acceptable formats");
	SUBREADputs("              \tinclude `GTF' and `SAF'. `GTF' by default. Please refer to the");
	SUBREADputs("              \tusers guide for SAF annotation format."); 
	SUBREADputs("    "); 
	SUBREADputs("    -t <input>\tSpecify the feature type. Only rows which have the matched"); 
	SUBREADputs("              \tmatched feature type in the provided GTF annotation file"); 
	SUBREADputs("              \twill be included for read counting. `exon' by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -g <input>\tSpecify the attribute type used to group features (eg. exons)");
	SUBREADputs("              \tinto meta-features (eg. genes), when GTF annotation is provided.");
	SUBREADputs("              \t`gene_id' by default. This attribute type is usually the gene");
	SUBREADputs("              \tidentifier. This argument is useful for the meta-feature level");
	SUBREADputs("              \tsummarization.");
	SUBREADputs("    "); 
	SUBREADputs("    -f        \tIf specified, read summarization will be performed at the "); 
	SUBREADputs("              \tfeature level (eg. exon level). Otherwise, it is performed at");
	SUBREADputs("              \tmeta-feature level (eg. gene level).");
	SUBREADputs("    "); 
	SUBREADputs("    -O        \tIf specified, reads (or fragments if -p is specified) will"); 
	SUBREADputs("              \tbe allowed to be assigned to more than one matched meta-"); 
	SUBREADputs("              \tfeature (or feature if -f is specified). "); 
	SUBREADputs("    "); 
	SUBREADputs("    -s <int>  \tIndicate if strand-specific read counting should be performed.");
	SUBREADputs("              \tIt has three possible values:  0 (unstranded), 1 (stranded) and");
	SUBREADputs("              \t2 (reversely stranded). 0 by default.");
	SUBREADputs("    "); 
	SUBREADputs("    -M        \tIf specified, multi-mapping reads/fragments will be counted (ie.");
	SUBREADputs("              \ta multi-mapping read will be counted up to N times if it has N");
	SUBREADputs("              \treported mapping locations). The program uses the `NH' tag to");
	SUBREADputs("              \tfind multi-mapping reads.");
	SUBREADputs("    "); 
	SUBREADputs("    -Q <int>  \tThe minimum mapping quality score a read must satisfy in order");
	SUBREADputs("              \tto be counted. For paired-end reads, at least one end should");
	SUBREADputs("              \tsatisfy this criteria. 0 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -T <int>  \tNumber of the threads. 1 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -R        \tOutput read counting result for each read/fragment. For each");
	SUBREADputs("              \tinput read file, read counting results for reads/fragments will");
	SUBREADputs("              \tbe saved to a tab-delimited file that contains four columns");
	SUBREADputs("              \tincluding read name, status(assigned or the reason if not");
	SUBREADputs("              \tassigned), name of target feature/meta-feature and number of");
	SUBREADputs("              \thits if the read/fragment is counted multiple times. Name of");
	SUBREADputs("              \tthe file is the same as name of the input read file except a");
	SUBREADputs("              \tsuffix `.featureCounts' is added.");
	SUBREADputs("    "); 
	SUBREADputs("    --primary \tIf specified, only primary alignments will be counted. Primary");
	SUBREADputs("              \tand secondary alignments are identified using bit 0x100 in the");
	SUBREADputs("              \tFlag field of SAM/BAM files. All primary alignments in a dataset");
	SUBREADputs("              \twill be counted no matter they are from multi-mapping reads or");
	SUBREADputs("              \tnot ('-M' is ignored). ");
	SUBREADputs("    "); 
	SUBREADputs("    --readExtension5 <int>      Reads are extended upstream by <int> bases from");
	SUBREADputs("              \ttheir 5' end."); 
	SUBREADputs("    "); 
	SUBREADputs("    --readExtension3 <int>      Reads are extended upstream by <int> bases from");
	SUBREADputs("              \ttheir 3' end."); 
	SUBREADputs("    "); 
	SUBREADputs("    --minReadOverlap <int>      Specify the minimum number of overlapped bases");
	SUBREADputs("              \trequired to assign a read to a feature. 1 by default. Negative");
	SUBREADputs("              \tvalues are permitted, indicating a gap being allowed between a");
	SUBREADputs("              \tread and a feature.");
	SUBREADputs("    "); 
	SUBREADputs("    --countSplitAlignmentsOnly  If specified, only split alignments (CIGAR");
	SUBREADputs("              \tstrings containing letter `N') will be counted. All the other");
	SUBREADputs("              \talignments will be ignored. An example of split alignments is");
	SUBREADputs("              \tthe exon-spanning reads in RNA-seq data.");
	SUBREADputs("    "); 
	SUBREADputs("    --read2pos <5:3>            The read is reduced to its 5' most base or 3'");
	SUBREADputs("              \tmost base. Read summarization is then performed based on the");
	SUBREADputs("              \tsingle base which the read is reduced to."); 
	SUBREADputs("    "); 
	SUBREADputs("    --ignoreDup                 If specified, reads that were marked as");
	SUBREADputs("              \tduplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM");
	SUBREADputs("              \tfile is used for identifying duplicate reads. In paired end");
	SUBREADputs("              \tdata, the entire read pair will be ignored if at least one end");
	SUBREADputs("              \tis found to be a duplicate read.");
	SUBREADputs("    "); 
	SUBREADputs("    Optional paired-end parameters:"); 
	SUBREADputs("    "); 
	SUBREADputs("    -p        \tIf specified, fragments (or templates) will be counted instead");
	SUBREADputs("              \tof reads. This option is only applicable for paired-end reads.");
	SUBREADputs("              \tThe two reads from the same fragment must be adjacent to each");
	SUBREADputs("              \tother in the provided SAM/BAM file.");
	SUBREADputs("    "); 
	SUBREADputs("    -P        \tIf specified, paired-end distance will be checked when assigning");
	SUBREADputs("              \tfragments to meta-features or features. This option is only");
	SUBREADputs("              \tapplicable when -p is specified. The distance thresholds should");
	SUBREADputs("              \tbe specified using -d and -D options."); 
	SUBREADputs("    "); 
	SUBREADputs("    -d <int>  \tMinimum fragment/template length, 50 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -D <int>  \tMaximum fragment/template length, 600 by default."); 
	SUBREADputs("    "); 
	SUBREADputs("    -B        \tIf specified, only fragments that have both ends "); 
	SUBREADputs("              \tsuccessfully aligned will be considered for summarization."); 
	SUBREADputs("              \tThis option is only applicable for paired-end reads."); 
	SUBREADputs("    "); 
	SUBREADputs("    -C        \tIf specified, the chimeric fragments (those fragments that "); 
	SUBREADputs("              \thave their two ends aligned to different chromosomes) will"); 
	SUBREADputs("              \tNOT be included for summarization. This option is only "); 
	SUBREADputs("              \tapplicable for paired-end read data."); 
	SUBREADputs("    "); 
	SUBREADputs("    -v        \tOutput version of the program.");
	SUBREADputs("    "); 


}



int readSummary_single_file(fc_thread_global_context_t * global_context, unsigned int * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter);

int readSummary(int argc,char *argv[]){

	/*
	   This function counts the number of reads falling into each exon region.
	   The order of exons in the output is the same as that of exons included in the annotation.
	   The annotation, if provided as a file, should be sorted by chromosome name.

	   Parameters passed from the featureCounts R function:
	0: "readSummary"
	1: ann
	2: files[i]
	3: fout
	4: as.numeric(isPairedEnd)
	5: min.distance
	6: max.distance
	7: as.numeric(tolower(file.type)=="sam")
	8: as.numeric(allowMultiOverlap)
	9: as.numeric(isGeneLevel)
	10: as.numeric(nthreads)
	11: as.numeric(isGTFannotation)
	12: as.numeric(isStrandChecked)
	13: as.numeric(isReadSummaryReported)
	14: as.numeric(isBothEndMapped)
	15: as.numeric(isChimericDisallowed)
	16: as.numeric(isPEDistChecked)
	17: nameFeatureTypeColumn 
	18: nameGeneIDColumn
	19: min.MappingQualityScore
	20: as.numeric(isMultiMappingAllowed)
	21: Annotation Chromosome Alias Name File. If the file is not specified, set this value to NULL or a zero-length string.
	22: Command line for CfeatureCounts header output; RfeatureCounts should set this value to NULL or a zero-length string or a space (' ').
	23: as.numeric(isInputFileResortNeeded)
	24: NOT IN USE: as.numeric(feature_block_size) # This parameter is no longer used. Give "14" for safe. 
	25: as.numeric(Five_End_Extension_Length)  # 5' end extension
	26: as.numeric(Three_End_Extension_Length)  # 3' end extension
	27: as.numeric(Minimum_Overlap_Between_Read_And_Feature) # 1 by default
	28: as.numeric(is_Split_Alignment_Only) # 0 by default
	29: as.numeric(reduce_5_3_ends_to_one) # 0= no reduction; 1= reduce to 5' end; 2= reduce to 3' end
	30: debug_command # This is for debug only; RfeatureCounts should pass a space (" ") to this parameter, disabling the debug command.
	31: as.numeric(is_duplicate_ignored) # 0 = INCLUDE DUPLICATE READS; 1 = IGNORE DUPLICATE READS (0x400 FLAG IS SET) ; "0" by default.
	 */

	int isStrandChecked, isCVersion, isChimericDisallowed, isPEDistChecked, minMappingQualityScore=0, isInputFileResortNeeded, feature_block_size = 20, reduce_5_3_ends_to_one;
	char **chr;
	long *start, *stop;
	int *geneid;

	char *nameFeatureTypeColumn, *nameGeneIDColumn,*debug_command;
	long nexons;


	long * anno_chr_head, * block_min_start, *block_max_end, *block_end_index;
	char ** anno_chrs, * anno_chr_2ch;
	long curchr, curpos;
	char * curchr_name;
	unsigned char * sorted_strand;
	curchr = 0;
	curpos = 0;
	curchr_name = "";

	int isPE, minPEDistance, maxPEDistance, isReadSummaryReport, isBothEndRequired, isMultiMappingAllowed, fiveEndExtension, threeEndExtension, minReadOverlap, isSplitAlignmentOnly, is_duplicate_ignored;

	int isSAM, isGTF, n_input_files=0;
	char *  alias_file_name = NULL, * cmd_rebuilt = NULL;

	int isMultiOverlapAllowed, isGeneLevel;

	isCVersion = ((argv[0][0]=='C')?1:0);

	isPE = atoi(argv[4]);
	minPEDistance = atoi(argv[5]);
	maxPEDistance = atoi(argv[6]);

	isSAM = atoi(argv[7]);
	isMultiOverlapAllowed = atoi(argv[8]);
	isGeneLevel = atoi(argv[9]);
	unsigned short thread_number;
	if(argc > 10)
		thread_number = atoi(argv[10]);
	else	thread_number = 4;
	if(argc > 11)
		isGTF = atoi(argv[11]);
	else	isGTF = 0;
	if(argc > 12)
		isStrandChecked = atoi(argv[12]);
	else	isStrandChecked = 0;
	if(argc > 13)
		isReadSummaryReport = atoi(argv[13]);
	else	isReadSummaryReport = 0;
	if(argc > 14)
		isBothEndRequired = atoi(argv[14]);
	else	isBothEndRequired = 0;
	if(argc > 15)
		isChimericDisallowed = atoi(argv[15]);
	else	isChimericDisallowed = 0;
	if(argc > 16)
		isPEDistChecked = atoi(argv[16]);
	else	isPEDistChecked = 0;
	if(argc > 17)
		nameFeatureTypeColumn = argv[17];
	else	nameFeatureTypeColumn = "exon";
	if(argc > 18)
		nameGeneIDColumn = argv[18];
	else	nameGeneIDColumn = "gene_id";
	if(argc > 19)
		minMappingQualityScore = atoi(argv[19]);
	else	minMappingQualityScore = 0;
	if(argc > 20)
		isMultiMappingAllowed = atoi(argv[20]);
	else	isMultiMappingAllowed = 1;
	if(argc > 21)
	{
		alias_file_name = argv[21];
		if(alias_file_name == NULL || alias_file_name[0]==' ' || alias_file_name[0]==0)
			alias_file_name = NULL;
	}
	else	alias_file_name = NULL;
	if(argc > 22)
	{
		cmd_rebuilt = argv[22];
		if(cmd_rebuilt == NULL || cmd_rebuilt[0]==' '||cmd_rebuilt[0]==0)
			cmd_rebuilt=NULL;
	}
	else	cmd_rebuilt = NULL;
	if(argc>23)
		isInputFileResortNeeded = atoi(argv[23]);
	else	isInputFileResortNeeded = 0;
	if(thread_number<1) thread_number=1;
	if(thread_number>16)thread_number=16;

	if(argc>25)
		fiveEndExtension = atoi(argv[25]);
	else 	fiveEndExtension = 0;

	if(argc>26)
		threeEndExtension = atoi(argv[26]);
	else	threeEndExtension = 0;

	if(argc>27)
		minReadOverlap = atoi(argv[27]);
	else	minReadOverlap = 1;
	
	if(argc>28)
		isSplitAlignmentOnly = atoi(argv[28]);
	else	isSplitAlignmentOnly = 0;

	if(argc>29)
		reduce_5_3_ends_to_one = atoi(argv[29]);	// 0 : no reduce; 1: reduce to 5' end; 2: reduce to 3' end.
	else	reduce_5_3_ends_to_one = 0;


	if(argc>30 && strlen(argv[30])>0 && argv[30][0]!=' ')
		debug_command = argv[30];
	else
		debug_command = " ";

	if(argc>31)
		is_duplicate_ignored = atoi(argv[31]);
	else
		is_duplicate_ignored = 0;


	unsigned int buffer_size = 1024*1024*6;


	fc_thread_global_context_t global_context;
	fc_thread_init_global_context(& global_context, buffer_size, thread_number, MAX_LINE_LENGTH, isPE, minPEDistance, maxPEDistance,isGeneLevel, isMultiOverlapAllowed, isStrandChecked, (char *)argv[3] , isReadSummaryReport, isBothEndRequired, isChimericDisallowed, isPEDistChecked, nameFeatureTypeColumn, nameGeneIDColumn, minMappingQualityScore,isMultiMappingAllowed, isSAM, alias_file_name, cmd_rebuilt, isInputFileResortNeeded, feature_block_size, isCVersion, fiveEndExtension, threeEndExtension , minReadOverlap, isSplitAlignmentOnly, reduce_5_3_ends_to_one, debug_command, is_duplicate_ignored);
	print_FC_configuration(&global_context, argv[1], argv[2], argv[3], global_context.is_SAM_file, isGTF, & n_input_files, isReadSummaryReport);


	//print_in_box(80,0,0,"IG=%d, IS=%d", isGeneLevel, isSplitAlignmentOnly);
	if(isSplitAlignmentOnly && ( isGeneLevel || !isMultiOverlapAllowed) )
	{
		print_in_box(80,0,0,"NOTICE --countSplitAlignmentsOnly is specified, but '-O' and '-f' are not");
		print_in_box(80,0,0,"       both specified. Please read the manual for details.");
		print_in_box(80,0,0,"");
	}

	// Loading the annotations.
	// Nothing is done if the annotation does not exist.
	fc_feature_info_t * loaded_features;
	print_in_box(84,0,0,"Load annotation file %s %c[0m...", argv[1], CHAR_ESC);
	nexons = load_feature_info(&global_context,argv[1], isGTF?FILE_TYPE_GTF:FILE_TYPE_RSUBREAD, &loaded_features);
	if(nexons<1){
		SUBREADprintf("Failed to open the annotation file %s, or its format is incorrect, or it contains no '%s' features.\n",argv[1], nameFeatureTypeColumn);
		return -1;
	}

	sort_feature_info(&global_context, nexons, loaded_features, &chr, &geneid, &start, &stop, &sorted_strand, &anno_chr_2ch, &anno_chrs, &anno_chr_head, & block_end_index, & block_min_start, & block_max_end);
	print_in_box(80,0,0,"   Meta-features : %d", global_context . gene_name_table -> numOfElements);
	print_in_box(80,0,0,"   Chromosomes : %d", global_context . exontable_nchrs);

	print_in_box(80,0,0,"");

	global_context.exontable_exons = nexons;
	unsigned int * nreads = (unsigned int *) calloc(nexons,sizeof(int));









	char * tmp_pntr = NULL;
	char * file_list_used = malloc(strlen(argv[2])+1);
	strcpy(file_list_used, argv[2]);
	char * next_fn = strtok_r(file_list_used,";", &tmp_pntr);
	unsigned int ** table_columns = calloc( n_input_files , sizeof(int *)), i_files=0;
	fc_read_counters * read_counters = calloc(n_input_files , sizeof(fc_read_counters)); 

	for(;;){
		int redoing, original_sorting = global_context.is_input_file_resort_needed, orininal_isPE = global_context.is_paired_end_data;
		if(next_fn==NULL || strlen(next_fn)<1) break;

		unsigned int * column_numbers = calloc(nexons, sizeof(unsigned int ));

		strcpy(global_context.input_file_name, next_fn);
		strcpy(global_context.raw_input_file_name, next_fn);
		global_context.redo=0;
		
		for(redoing = 0; redoing < 1 + !original_sorting; redoing++)
		{
			fc_read_counters * my_read_counter = &(read_counters[i_files]);
			memset(my_read_counter, 0, sizeof(fc_read_counters));

			int ret_int = readSummary_single_file(& global_context, column_numbers, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start, block_max_end, my_read_counter);
			if(ret_int!=0 || (global_context.redo && redoing)){
				table_columns[i_files] = NULL;
				free(column_numbers);
				break;
			}
			else table_columns[i_files] = column_numbers;

			if(redoing || !global_context.redo) break;
			
			global_context.is_input_file_resort_needed = 1;
			memset(column_numbers, 0, nexons * sizeof(unsigned int ));
		}
		global_context.is_SAM_file = isSAM;
		global_context.is_input_file_resort_needed = original_sorting;
		global_context.is_paired_end_data = orininal_isPE;

		i_files++;
		next_fn = strtok_r(NULL, ";", &tmp_pntr);
	}

	free(file_list_used);

	if(isGeneLevel)
		fc_write_final_gene_results(&global_context, geneid, chr, start, stop, sorted_strand, argv[3], nexons,  table_columns, argv[2], n_input_files , loaded_features, isCVersion);
	else
		fc_write_final_results(&global_context, argv[3], nexons, table_columns, argv[2], n_input_files ,loaded_features, isCVersion);

	fc_write_final_counts(&global_context, argv[3], n_input_files, argv[2], table_columns, read_counters, isCVersion);

	for(i_files=0; i_files<n_input_files; i_files++)
		if(table_columns[i_files]) free(table_columns[i_files]);
	free(table_columns);


	print_FC_results(&global_context);
	KeyValuePair * cursor;
	int bucket;
	for(bucket=0; bucket < global_context.exontable_chro_table  -> numOfBuckets; bucket++)
	{
		cursor = global_context.exontable_chro_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			fc_chromosome_index_info * del_chro_info = cursor->value;
			free(del_chro_info->reverse_table_start_index);
			//free(del_chro_info->reverse_table_end_index);
			free((void *)cursor -> key);
			free(del_chro_info);
			cursor = cursor->next;
		}
	}

	if(global_context.SAM_output_fp) fclose(global_context. SAM_output_fp);
	HashTableDestroy(global_context.gene_name_table);
	free(global_context.gene_name_array);

	HashTableDestroy(global_context.exontable_chro_table);
	if(global_context.annot_chro_name_alias_table)
		HashTableDestroy(global_context.annot_chro_name_alias_table);

	free(global_context.unistr_buffer_space);
	free(loaded_features);
	free(geneid);
	free(chr);
	free(start);
	free(sorted_strand);
	free(anno_chr_2ch);
	free(anno_chrs);
	free(anno_chr_head);
	free(block_min_start);
	free(block_max_end);
	free(block_end_index);
	free(stop);
	free(nreads);


	return 0;
}













int readSummary_single_file(fc_thread_global_context_t * global_context, unsigned int * column_numbers, int nexons,  int * geneid, char ** chr, long * start, long * stop, unsigned char * sorted_strand, char * anno_chr_2ch, char ** anno_chrs, long * anno_chr_head, long * block_end_index, long * block_min_start , long * block_max_end, fc_read_counters * my_read_counter)
{
	FILE *fp_in = NULL;
	int read_length = 0;
	int is_first_read_PE=0;
	char * line = (char*)calloc(MAX_LINE_LENGTH, 1);
	char * file_str = "";

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		int file_probe = is_certainly_bam_file(global_context->input_file_name, &is_first_read_PE);

		// a Singel-end SAM/BAM file cannot be assigned as a PE SAM/BAM file;
		// but a PE SAM/BAM file may be assigned as a SE file if the user wishes to do so.
		if(is_first_read_PE==0)
				global_context -> is_paired_end_data = 0;

		if(file_probe == 1){
			global_context->is_SAM_file = 0;
		}
		else if(file_probe == 0) global_context->is_SAM_file = 1;

		global_context -> start_time = miltime();

		file_str = "SAM";
		if(file_probe == 1) file_str = "BAM" ;
		if(file_probe == -1) file_str = "Unknown";

		if(!global_context->redo)
		{
			print_in_box(80,0,0,"Process %s file %s...", file_str, global_context->input_file_name);
			if(is_first_read_PE)
				print_in_box(80,0,0,"   Paired-end reads are included.");
			else
				print_in_box(80,0,0,"   Single-end reads are included.");
		}
		
	}

	int isInputFileResortNeeded = global_context->is_input_file_resort_needed;

	if(strcmp( global_context->input_file_name,"STDIN")!=0)
	{
		FILE * exist_fp = f_subr_open( global_context->input_file_name,"r");
		if(!exist_fp)
		{
			print_in_box(80,0,0,"Failed to open file %s",  global_context->input_file_name);
			print_in_box(80,0,0,"No counts were generated for this file.");
			print_in_box(80,0,0,"");
			return -1;
		}
		fclose(exist_fp);
	}

	if(strcmp(global_context->input_file_name,"STDIN")!=0)
		if(warning_file_type(global_context->input_file_name, global_context->is_SAM_file?FILE_TYPE_SAM:FILE_TYPE_BAM))
			global_context->is_unpaired_warning_shown=1;
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		if(resort_input_file( global_context)) return -1;
	int isSAM = global_context->is_SAM_file;
	// Open the SAM/BAM file
	// Nothing is done if the file does not exist.

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)==0)
		fp_in = stdin;
	else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#else
		fp_in = f_subr_open(global_context->input_file_name,"r");
	#endif


	// begin to load-in the data.
	if(!global_context->redo)
	{
		if( global_context->is_paired_end_data)
		{
			print_in_box(80,0,0,"   Assign fragments (read pairs) to features...");
//				print_in_box(80,0,0,"   Each fragment is counted no more than once.");
		}
		else
			print_in_box(80,0,0,"   Assign reads to features...");
	}



	fc_thread_start_threads(global_context, nexons, geneid, chr, start, stop, sorted_strand, anno_chr_2ch, anno_chrs, anno_chr_head, block_end_index, block_min_start , block_max_end, read_length);

	int buffer_pairs = global_context -> thread_number>1?512:1;
	int isPE = global_context->is_paired_end_data;
	char * preload_line = malloc(sizeof(char) * (2+MAX_LINE_LENGTH)*(isPE?2:1)*buffer_pairs);
	int preload_line_ptr;
	int current_thread_id = 0;
	fc_thread_thread_context_t * one_thread_context = global_context->thread_contexts;

	SamBam_Reference_Info * sb_header_tab = NULL;
	
	unsigned long long int chunk_id = 0;
	int binary_remainder = 0, binary_read_ptr = 0;
	char * chunk_in_buff = malloc(70000);
	char * binary_in_buff = malloc(80000 * 2);

	if(!isSAM)
	{
		int remainder_read_data_len = 0;

		PBum_load_header(fp_in, &sb_header_tab, binary_in_buff,  & remainder_read_data_len);
		//printf("RMD=%d\n", remainder_read_data_len);

		if(remainder_read_data_len)
		{
			binary_remainder = remainder_read_data_len;
		}
		global_context->sambam_chro_table = sb_header_tab;
	}

	while (1){
		int pair_no;
		int is_second_read;
		int fresh_read_no = 0;
		preload_line[0] = 0;
		preload_line_ptr = 0;

		char * ret = NULL;
		
		// one-thread BAM is not supported.
		if( isSAM && global_context->thread_number==1)
		{
			int is_second_read;

			for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
			{
				char * lbuf = is_second_read?one_thread_context -> line_buffer2:one_thread_context -> line_buffer1;
				while(1)
				{
					ret = fgets(lbuf, MAX_LINE_LENGTH, fp_in);
					if(global_context -> redo) ret = NULL;
					if(!ret) break;
					if(lbuf[0] == '@') 
					{
						int retlen = strlen(ret);
						if(ret[retlen-1]!='\n')
						{
							while(1){
								int nch = getc(fp_in);
								if(nch == EOF || nch == '\n') break;
							}
						}
					}
					else break;
				}

				if(!ret) break;
				if(read_length < 1)
				{
					int tab_no = 0;
					int read_len_tmp=0, read_cursor;
					int curr_line_len = strlen(lbuf);
					for(read_cursor=0; read_cursor<curr_line_len; read_cursor++)
					{
						if(lbuf[read_cursor] == '\t')
							tab_no++;
						else
						{
							if(tab_no == 9)	// SEQ
								read_len_tmp++;
						}
					}
					read_length = read_len_tmp;
					global_context->read_length = read_length;
				}
			}

			//printf("RRR=%d\n",ret);
			if(!ret) break;
			
			//one_thread_context -> current_read_length1 = global_context->read_length;
			//one_thread_context -> current_read_length2 = global_context->read_length;

			global_context->all_reads ++;
			process_line_buffer(global_context, one_thread_context);
		}
		else if(!isSAM)
		{
			int no_of_reads = 0;
			unsigned int real_len = 0;
			// most of the data must have been given out before this step.

			int cdata_size = 0;

			if(binary_remainder > 70000)
				SUBREADprintf("SOMETHING IS WRONG!\n");

			if(global_context -> redo)
				cdata_size = -1;
			else{
				if(binary_remainder<10000)
					cdata_size = PBam_get_next_zchunk(fp_in, chunk_in_buff, 65537, & real_len);
			}

			if(cdata_size>0 || binary_remainder>0)
			{
				int x1;


				if(binary_read_ptr>0)
				{
					for(x1=0; x1< binary_remainder; x1++)
						binary_in_buff[x1] = binary_in_buff [x1 + binary_read_ptr];
					binary_read_ptr = 0;
				}

				//fprintf(stderr,"NBN=%d, OBN=%d\n", cdata_size , binary_remainder);
				if(cdata_size>0)
				{
					int new_binary_bytes = SamBam_unzip(binary_in_buff + binary_remainder , chunk_in_buff , cdata_size);
					if(new_binary_bytes>=0)
						binary_remainder += new_binary_bytes;
					else	SUBREADprintf("ERROR: BAM GZIP FORMAT ERROR.\n");
				//	fprintf(stderr,"BBN=%d\n", new_binary_bytes);
				}

				while(binary_remainder>4)
				{
					unsigned int binary_read_len = 0;
					memcpy(& binary_read_len , binary_in_buff + binary_read_ptr , 4);
					//printf("RLEN=%d; PTR=%d; RMD=%d\n", binary_read_len , binary_read_ptr, binary_remainder);
					if(binary_read_len > 10000)
					{
						binary_remainder = -1;
						//SUBREADprintf("FATAL ERROR: BAM RECORD SIZE = %u ; PTR=%d ; REM=%d.\n", binary_read_len, binary_read_ptr , binary_remainder);
						print_in_box(80,0,0,"   A format error was detected in this BAM file.");
						print_in_box(80,0,0,"   The remaining part in the file is skipped.");
						print_in_box(80,0,0,"   Please check the file format using samtools.");
						print_in_box(80,0,0,"");
						break;
					}
					// if the program runs on PE mode, no_of_reads must be even.

					if(isPE)
					{
						if(4 + binary_read_len + 4 < binary_remainder)
						{
							int binary_read2_len=0;
							memcpy(&binary_read2_len , binary_in_buff + binary_read_ptr + 4 + binary_read_len, 4);
							if(4 + binary_read_len + 4 + binary_read2_len <= binary_remainder)
							{
								no_of_reads +=2;
								binary_read_ptr += 4 + binary_read_len + 4 + binary_read2_len;
								binary_remainder -= 4 + binary_read_len + 4 + binary_read2_len;
							}
							else break;
						}
						else break;
					}
					else
					{
						if(binary_read_len + 4<= binary_remainder)
						{
							no_of_reads ++;
							binary_read_ptr  += 4 + binary_read_len;
							binary_remainder -= 4 + binary_read_len;
						}
						else break;
					}
				}
			}

			if(binary_remainder <0)break;

			if(no_of_reads>0) 
			{
				while(1)
				{
					int is_finished = 0;

					fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;

					pthread_spin_lock(&thread_context->input_buffer_lock);

					// the number of bytes can be utilised given the two_chunk_len.
					int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder;
					int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr;
					
					if(thread_context->input_buffer_remainder > global_context->input_buffer_max_size)
					{
						SUBREADprintf("RMD=%d\n", thread_context->input_buffer_remainder );
						assert(0);
					}

					if(tail_bytes < binary_read_ptr + 4)
						empty_bytes -= tail_bytes;

					// copy the new buffer to thread buffer.
					// format: read_number=n, read_chunk1, read_chunk2, ..., read_chunk_n
					if(empty_bytes >= binary_read_ptr + 4)
					{

						if(tail_bytes < binary_read_ptr + 4)
						{
							if(tail_bytes>=4)
								memset(thread_context->input_buffer + thread_context->input_buffer_write_ptr, 0, 4);
							thread_context->input_buffer_write_ptr = 0;
							thread_context->input_buffer_remainder += tail_bytes;
						}

						memcpy( thread_context->input_buffer + thread_context->input_buffer_write_ptr, & binary_read_ptr, 4);
						memcpy( thread_context->input_buffer + thread_context->input_buffer_write_ptr + 4, binary_in_buff , binary_read_ptr);
						thread_context->input_buffer_write_ptr += 4 + binary_read_ptr;
						thread_context->input_buffer_remainder += 4 + binary_read_ptr;
						is_finished = 1;
					}

					pthread_spin_unlock(&thread_context->input_buffer_lock);
					current_thread_id++;
					if(current_thread_id >= global_context->thread_number) current_thread_id = 0;

					if(is_finished) break;
					else usleep(tick_time);

				}

				chunk_id++;
			}
			else if(cdata_size<0)
				break;
	
		}
		else
		{
			for(pair_no=0; pair_no < buffer_pairs; pair_no++)
			{
				for(is_second_read=0;is_second_read<(isPE?2:1);is_second_read++)
				{
					while(1)
					{
						ret = fgets(preload_line+preload_line_ptr, MAX_LINE_LENGTH, fp_in);
						if(global_context -> redo ) ret = NULL;
						if(!ret) break;
						if(preload_line[preload_line_ptr] == '@'){
							int retlen = strlen(ret);
							if(ret[retlen-1]!='\n')
							{
								while(1){
									int nch = getc(fp_in);
									if(nch == EOF || nch == '\n') break;
								}
							}
						}else break;
					}

					if(!ret) break;

					int curr_line_len = strlen(preload_line+preload_line_ptr);
					if(curr_line_len >= MAX_LINE_LENGTH || preload_line[preload_line_ptr + curr_line_len-1]!='\n')
					{
						print_in_box(80,0,0,"ERROR: the lines are too long. Please check the input format!!\n");
						ret = NULL;
						preload_line_ptr = 0;
						break;
					}
					preload_line_ptr += curr_line_len;

					fresh_read_no++;
				}
				if(!ret) break;
				else
					global_context->all_reads ++;
			}

			int line_length = preload_line_ptr;
			if(line_length >= global_context->input_buffer_max_size-1)
			{
				SUBREADprintf("ERROR: the lines are too long. Please check the input format!!\n");
				break;
			}
			if(isPE && (fresh_read_no%2>0))
			{
				// Safegarding -- it should not happen if the SAM file has a correct format.
				line_length = 0;
				if( (!global_context -> redo)){
					print_in_box(80,0,0,"   There are odd number of reads in the paired-end data.");
					print_in_box(80,0,0,"   Please make sure that the format is correct.");
				}
			}

			//printf("FRR=%d\n%s\n", fresh_read_no, preload_line);

			if(line_length > 0)
			{
				while(1)
				{
					int is_finished = 0;
					fc_thread_thread_context_t * thread_context = global_context->thread_contexts+current_thread_id;
					//printf("WRT_THR_IBUF_REM [%d]=%d\n", current_thread_id , thread_context->input_buffer_remainder);

					pthread_spin_lock(&thread_context->input_buffer_lock);
					unsigned int empty_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_remainder; 
					if(empty_bytes > line_length)
					{
						unsigned int tail_bytes = global_context->input_buffer_max_size -  thread_context->input_buffer_write_ptr; 
						unsigned int write_p1_len = (tail_bytes > line_length)?line_length:tail_bytes;
						unsigned int write_p2_len = (tail_bytes > line_length)?0:(line_length - tail_bytes);
						memcpy(thread_context->input_buffer + thread_context->input_buffer_write_ptr, preload_line, write_p1_len);
						if(write_p2_len)
						{
							memcpy(thread_context->input_buffer, preload_line + write_p1_len, write_p2_len);
							thread_context->input_buffer_write_ptr = write_p2_len;
						}
						else	thread_context->input_buffer_write_ptr += write_p1_len;
						if(thread_context->input_buffer_write_ptr == global_context->input_buffer_max_size) 
							thread_context->input_buffer_write_ptr=0;


						thread_context->input_buffer_remainder += line_length;
						//printf("WRT_THR_IBUF_REM [%d] + %d =%d\n", current_thread_id, line_length , thread_context->input_buffer_remainder);
						is_finished = 1;
					}

					pthread_spin_unlock(&thread_context->input_buffer_lock);

					current_thread_id++;
					if(current_thread_id >= global_context->thread_number) current_thread_id = 0;

					if(is_finished) break;
					else usleep(tick_time);
				}
			}
			if(!ret) break;
		}
	}

	free(chunk_in_buff);
	free(binary_in_buff);
	free(preload_line);
	global_context->is_all_finished = 1;

	if(global_context->thread_number > 1 || !isSAM)
		fc_thread_wait_threads(global_context);

	unsigned long long int nreads_mapped_to_exon = 0;


	if(!global_context->redo)
		fc_thread_merge_results(global_context, column_numbers , &nreads_mapped_to_exon, my_read_counter);

	fc_thread_destroy_thread_context(global_context);

	//global_context .read_counters.assigned_reads = nreads_mapped_to_exon;

	#ifdef MAKE_STANDALONE
	if(strcmp("STDIN",global_context->input_file_name)!=0)
	#endif
		fclose(fp_in);

	if(sb_header_tab) free(sb_header_tab);
	if(strcmp(global_context->input_file_name,"STDIN")!=0 && isInputFileResortNeeded)
		unlink(global_context->input_file_name);
	free(line);
	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int feature_count_main(int argc, char ** argv)
#endif
{
	char * Rargv[32];
	char annot_name[300];
	char * out_name = malloc(300);
	char * alias_file_name = malloc(300);
	int cmd_rebuilt_size = 200;
	char * cmd_rebuilt = malloc(cmd_rebuilt_size);
	char nameFeatureTypeColumn[66];
	char nameGeneIDColumn[66];
	int min_qual_score = 0;
	int min_dist = 50;
	int max_dist = 600;
	char debug_command[10];
	char min_dist_str[11];
	char max_dist_str[11];
	char min_qual_score_str[11];
	char feature_block_size_str[11];
	char Strand_Sensitive_Str[11];
	char * very_long_file_names;
	int is_Input_Need_Reorder = 0;
	int is_PE = 0;
	int is_SAM = 1;
	int is_GeneLevel = 1;
	int is_Overlap = 0;
	int is_Both_End_Mapped = 0;
	int feature_block_size = 14;
	int Strand_Sensitive_Mode = 0;
	int is_ReadSummary_Report = 0;
	int is_Chimeric_Disallowed = 0;
	int is_PE_Dist_Checked = 0;
	int is_Multi_Mapping_Allowed = 0;
	int is_Split_Alignment_Only = 0;
	int is_duplicate_ignored = 0;
	int reduce_5_3_ends_to_one = 0;
	int threads = 1;
	int isGTF = 1;
	char nthread_str[4];
	int option_index = 0;
	int c;
	int very_long_file_names_size = 200;
	int fiveEndExtension = 0, threeEndExtension = 0, minReadOverlap = 1;
	char strFiveEndExtension[11], strThreeEndExtension[11], strMinReadOverlap[11];
	very_long_file_names = malloc(very_long_file_names_size);
	very_long_file_names [0] = 0;

	alias_file_name[0]=0;
	debug_command[0] = 0;

	strcpy(nameFeatureTypeColumn,"exon");
	strcpy(nameGeneIDColumn,"gene_id");
	annot_name[0]=0;out_name[0]=0;
	

	cmd_rebuilt[0]=0;
	for(c = 0; c<argc;c++)
	{
		if(strlen(cmd_rebuilt) + 300 > cmd_rebuilt_size)
		{
			cmd_rebuilt_size*=2;
			cmd_rebuilt = realloc(cmd_rebuilt, cmd_rebuilt_size);
		}
		sprintf(cmd_rebuilt+strlen(cmd_rebuilt), "\"%s\" ", argv[c]);
	}

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "A:g:t:T:o:a:d:D:L:Q:pbF:fs:SCBPMORv?", long_options, &option_index)) != -1)
		switch(c)
		{
			case 'S':
				SUBREADprintf("The '-S' option has been deprecated.\n FeatureCounts will automatically examine the read order.\n");
				//is_Input_Need_Reorder = 1;
				break;
			case 'A':
				strcpy(alias_file_name, optarg);
				break;
			case 'M':
				if(0 == is_Multi_Mapping_Allowed)
					is_Multi_Mapping_Allowed = ALLOW_ALL_MULTI_MAPPING;
				break;
			case 'v':
				core_version_number("featureCounts");
				return 0;
			case 'Q':
				min_qual_score = atoi(optarg);
				break;
			case 't':
				strcpy(nameFeatureTypeColumn, optarg);
				break;
			case 'g':
				while((*optarg) == ' ') optarg++;
				strcpy(nameGeneIDColumn, optarg);
				break;
			case 'T':
				threads = atoi(optarg);
				break;
			case 'd':
				min_dist = atoi(optarg);
				break;
			case 'D':
				max_dist = atoi(optarg);
				break;
			case 'p':
				is_PE = 1;
				break;
			case 'b':
				SUBREADprintf("The '-b' option has been deprecated.\n FeatureCounts will automatically examine the file format.\n");
				is_SAM = 0;
				break;
			case 'C':
				is_Chimeric_Disallowed = 1;
				break;
			case 'P':
				is_PE_Dist_Checked = 1;
				break;
			case 'B':
				is_Both_End_Mapped = 1;
				break;
			case 'f':
				is_GeneLevel = 0;
				break;
			case 'F':
				isGTF = 1;
				if(strcmp("SAF", optarg)==0) isGTF=0;
				else if(strcmp("GTF", optarg)==0) isGTF=1;
				else SUBREADprintf("\nWarning: Unknown annotation format: %s. GTF format is used.\n\n", optarg); 
				break;
			case 'O':
				is_Overlap = 1;
				break;
			case 'R':
				is_ReadSummary_Report = 1;
				break;
			case 's':
				Strand_Sensitive_Mode = atoi(optarg);
				break;
//			case 'i':
//				term_strncpy(sam_name, optarg,299);
//				break;
			case 'o':
				term_strncpy(out_name, optarg,299);
				break;
			case 'a':
				term_strncpy(annot_name, optarg,299);
				break;
			case 'L':
				feature_block_size = atoi(optarg);
				break;
			case 0 :	// long options

				if(strcmp("primary", long_options[option_index].name)==0)
				{
					is_Multi_Mapping_Allowed = ALLOW_PRIMARY_MAPPING;
				}

				if(strcmp("readExtension5", long_options[option_index].name)==0)
				{
					fiveEndExtension = atoi(optarg);
					fiveEndExtension = max(0, fiveEndExtension);
				}

				if(strcmp("readExtension3", long_options[option_index].name)==0)
				{
					threeEndExtension = atoi(optarg);
					threeEndExtension = max(0, threeEndExtension);
				}

				if(strcmp("minReadOverlap", long_options[option_index].name)==0)
				{
					minReadOverlap = atoi(optarg);
				}

				if(strcmp("debugCommand", long_options[option_index].name)==0)
				{
					strcpy(debug_command, optarg);
				}


				if(strcmp("ignoreDup", long_options[option_index].name)==0)
				{
					is_duplicate_ignored = 1 ;
				}

				if(strcmp("read2pos", long_options[option_index].name)==0)
				{
					if(optarg[0]=='3')
						reduce_5_3_ends_to_one = REDUCE_TO_3_PRIME_END;
					else if(optarg[0]=='5')
						reduce_5_3_ends_to_one = REDUCE_TO_5_PRIME_END;
						
				}				

				if(strcmp("countSplitAlignmentsOnly", long_options[option_index].name)==0)
				{
					is_Split_Alignment_Only = 1;
				}

				break;
			case '?':
			default :
				print_usage();
				return -1;
				break;
		}


	if(minReadOverlap<1)
	{
		fiveEndExtension = - minReadOverlap + 1;
		threeEndExtension =  - minReadOverlap + 1;
		minReadOverlap = 1;
	}

	if(out_name[0]==0 || annot_name[0]==0||argc == optind)
	{
		print_usage();
		return -1;
	}

	for(; optind < argc; optind++)
	{
		int curr_strlen = strlen(very_long_file_names);
		if( very_long_file_names_size - curr_strlen <300)
		{
			very_long_file_names_size *=2;
			//printf("CL=%d ; NS=%d\n", curr_strlen , very_long_file_names_size);
			very_long_file_names=realloc(very_long_file_names , very_long_file_names_size);
		}

		strcat(very_long_file_names, argv[optind]);
		strcat(very_long_file_names, ";");
	}

	very_long_file_names[strlen(very_long_file_names)-1]=0;

	sprintf(strFiveEndExtension, "%d", fiveEndExtension);
	sprintf(strThreeEndExtension, "%d", threeEndExtension);
	sprintf(strMinReadOverlap, "%d", minReadOverlap);
	sprintf(nthread_str,"%d", threads);
	sprintf(min_dist_str,"%d",min_dist);
	sprintf(max_dist_str,"%d",max_dist);
	sprintf(min_qual_score_str,"%d", min_qual_score);
	sprintf(feature_block_size_str,"%d", feature_block_size);
	sprintf(Strand_Sensitive_Str,"%d", Strand_Sensitive_Mode);
	Rargv[0] = "CreadSummary";
	Rargv[1] = annot_name;
	Rargv[2] = very_long_file_names;
	Rargv[3] = out_name;
	Rargv[4] = is_PE?"1":"0";
	Rargv[5] = min_dist_str;
	Rargv[6] = max_dist_str;
	Rargv[7] = is_SAM?"1":"0";
	Rargv[8] = is_Overlap?"1":"0";
	Rargv[9] = is_GeneLevel?"1":"0";
	Rargv[10] = nthread_str;
	Rargv[11] = isGTF?"1":"0";
	Rargv[12] = Strand_Sensitive_Str;
	Rargv[13] = is_ReadSummary_Report?"1":"0";
	Rargv[14] = is_Both_End_Mapped?"1":"0";
	Rargv[15] = is_Chimeric_Disallowed?"1":"0";
	Rargv[16] = is_PE_Dist_Checked?"1":"0";
	Rargv[17] = nameFeatureTypeColumn;
	Rargv[18] = nameGeneIDColumn;
	Rargv[19] = min_qual_score_str;
	Rargv[20] = is_Multi_Mapping_Allowed == ALLOW_PRIMARY_MAPPING?"2":(is_Multi_Mapping_Allowed == ALLOW_ALL_MULTI_MAPPING?"1":"0");
	Rargv[21] = alias_file_name;
	Rargv[22] = cmd_rebuilt;
	Rargv[23] = is_Input_Need_Reorder?"1":"0";
	Rargv[24] = feature_block_size_str;
	Rargv[25] = strFiveEndExtension;
	Rargv[26] = strThreeEndExtension;
	Rargv[27] = strMinReadOverlap;
	Rargv[28] = is_Split_Alignment_Only?"1":"0";
	Rargv[29] = (reduce_5_3_ends_to_one == 0?"0":(reduce_5_3_ends_to_one==REDUCE_TO_3_PRIME_END?"3":"5"));
	Rargv[30] = debug_command;
	Rargv[31] = is_duplicate_ignored?"1":"0";
	readSummary(32, Rargv);

	free(very_long_file_names);
	free(out_name);
	free(alias_file_name);
	free(cmd_rebuilt);

	return 0;

}


