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
  
  
#ifndef _SUBREAD_H_
#define _SUBREAD_H_

#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>

#ifndef MAKE_STANDALONE
#ifndef RUNNING_ENV
#include <R.h>
#endif
#endif

#include "hashtable.h" 

#define SAM_FLAG_PAIRED_TASK	0x01
#define SAM_FLAG_FIRST_READ_IN_PAIR 0x40
#define SAM_FLAG_SECOND_READ_IN_PAIR 0x80
#define SAM_FLAG_MATE_UNMATCHED 0x08
#define SAM_FLAG_MATCHED_IN_PAIR 0x02
#define SAM_FLAG_REVERSE_STRAND_MATCHED 0x10
#define SAM_FLAG_MATE_REVERSE_STRAND_MATCHED 0x20
#define SAM_FLAG_SECONDARY_MAPPING 0x100
#define SAM_FLAG_DUPLICATE 0x400
#define SAM_FLAG_UNMAPPED 0x04

#define FUSION_BREAK_POINT	2
#define FUSION_JUNCTION		1
#define SPLICING_JUNCTION	0

#define RUN_ALIGN 		0
#define RUN_FINAL 		1



#define MAX_PIECE_JUNCTION_READ 7
#define MAX_READ_LENGTH 1210
#define MAX_READ_NAME_LEN 100 
#define MAX_CHROMOSOME_NAME_LEN 100 
#define MAX_FILE_NAME_LENGTH 300

#define EXON_LONG_READ_LENGTH 120
#define EXON_MAX_CIGAR_LEN 48
#define MAX_INDEL_SECTIONS 7
//#define XBIG_MARGIN_RECORD_SIZE 24 
#define MAX_INSERTION_LENGTH 200
//#define BASE_BLOCK_LENGTH 15000000
//#define NEED_SUBREAD_STATISTIC


#define IS_MIN_POS_NEGATIVE_STRAND 4
#define IS_MAX_POS_NEGATIVE_STRAND 12 
#define IS_PAIRED_HINTED 16
#define IS_R1_CLOSE_TO_5 1
#define IS_REVERSED_HALVES 2
#define	IS_PROCESSED_READ 32
#define	IS_PROCESSED_READ_R2 64
#define IS_PAIRED_MATCH 128
#define IS_NEGATIVE_STRAND_R1 256 
#define IS_NEGATIVE_STRAND_R2 512 
#define IS_FUSION 1024 
#define IS_NEGATIVE_STRAND 2048
#define IS_RECOVERED_JUNCTION_READ 4096
#define IS_FINALISED_PROCESSING 8192
#define IS_RECOVERED_JUNCTION_READ_STEP4 (8192*2)
#define	IS_BREAKEVEN_READ (8192*4)
#define IS_R1R2_EQUAL_LEN 1024

#if defined(MACOS) || defined(FREEBSD)
typedef pthread_mutex_t subread_lock_t;
#define pthread_spinlock_t pthread_mutex_t
#define pthread_spin_lock pthread_mutex_lock
#define pthread_spin_unlock pthread_mutex_unlock
#define pthread_spin_init(a, b) pthread_mutex_init(a, NULL)
#define pthread_spin_destroy(a) pthread_mutex_destroy(a) 
#define strnlen(a,l) strlen(a)
#else
typedef pthread_spinlock_t subread_lock_t;
#endif

#ifdef WINDOWS
#define ftello ftello64
#define fseeko fseeko64

#endif

#if defined(MAKE_STANDALONE) || defined(RUNNING_ENV)
#define SUBREADprintf(...) fprintf(stderr, __VA_ARGS__)
#define SUBREADputs(x) fprintf(stderr, "%s\n", x)
#define SUBREADputchar(x) fputc(x, stderr)
#define SUBREADfflush(x) fflush(x)
#define CORE_SOFT_BR_CHAR '\n'
#else
#define SUBREADprintf Rprintf
#define SUBREADputs(x) Rprintf("%s\n",(x))
#define SUBREADputchar(X) Rprintf("%c",(X)) 
#define SUBREADfflush(X) 
#define CORE_SOFT_BR_CHAR '\n'

#endif

#ifndef NONONO_DONOTDEF

#define QUALITY_KILL	198
#define QUALITY_KILL_SUBREAD	160
#define MAX_QUALITY_TO_CALL_JUNCTION 2195
#define MAX_QUALITY_TO_EXPLORER_JUNCTION 209

#else

#define TEST_TARGET ""

#endif

#define SNP_CALLING_ONLY_HIGHQUAL 1

#define MESSAGE_OUT_OF_MEMORY "Out of memory. If you are using Rsubread in R, please save your working environment and restart R. \n"
#define fatal_memory_size(a) SUBREADputs(MESSAGE_OUT_OF_MEMORY);

//#define QUALITY_KILL	175
//#define QUALITY_KILL_SUBREAD	150



typedef unsigned int gehash_key_t;
typedef unsigned int gehash_data_t;
//typedef float gene_quality_score_t;
typedef int gene_quality_score_t;
//typedef unsigned char gene_vote_number_t;
typedef short gene_vote_number_t;


#define XOFFSET_TABLE_SIZE 250000

#define ANCHORS_NUMBER 259

#define BEXT_RESULT_LIMIT 16

#define SEARCH_BACK 0
#define SEARCH_FRONT 1

//#define LARGE_GENE_VOTE_TABLE
#ifdef LARGE_GENE_VOTE_TABLE
#warning "Using LARGE_GENE_VOTE_TABLE"
#define GENE_VOTE_SPACE 32 
#define GENE_VOTE_TABLE_SIZE 109 
#else
#define GENE_VOTE_SPACE 8
#define GENE_VOTE_TABLE_SIZE 61
#endif

#define MAX_ANNOTATION_EXONS 30000 
#define MAX_EXONS_PER_GENE 400 
#define MAX_EXON_CONNECTIONS 10

#define MAX_GENE_NAME_LEN 12
#define MAX_INDEL_TOLERANCE 7

#define SUBINDEX_VER0 100
#define SUBINDEX_VER1 200
#define SUBINDEX_VER2 201

#define SUBREAD_INDEX_OPTION_INDEX_GAP 0x0101


#define CHAR_ESC 27

//#define base2int(c) ((c)=='A'?0:((c)=='T'?3:((c)=='C'?2:1)))
#define base2int(c) ((c)<'G'?((c)=='A'?0:2):((c)=='G'?1:3))


                     // A  B  C  D  E  F  G
//#define base2int(c) ("\x0\x0\x2\x0\x0\x0\x1\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3\x3"[(c)-'A'])
//#define int2base(c) ((c)==1?'G':((c)==0?'A':((c)==2?'C':'T')))
//#define int2base(c) ("AGCT"[(c)]) 
#define int2base(c) (1413695297 >> (8*(c))&0xff);
#define color2int(c) ((c) - '0')
#define int2color(c) ("0123"[(c)])
#define remove_backslash(str) { int xxxa=0; while(str[xxxa]){ if(str[xxxa]=='/'){str[xxxa]='\0'; break;} xxxa++;} }

/*
#define get_base_error_prob64(a) ((a) < '@'-1?1:pow(10., -0.1*((a)-'@')))
#define get_base_error_prob33(a) ((a) < '!'-1?1:pow(10., -0.1*((a)-'!'))) 

*/
#define SUBREAD_malloc(a) malloc(a)

#define FASTQ_PHRED33 1
#define FASTQ_PHRED64 0

#define IS_DEBUG 0


typedef struct {
  char gene_name [MAX_GENE_NAME_LEN]; 
  // The chromosome name is not stored in this data structure
  // All coordinates are translated into the linear location in the entire referenced genome, usually 0 ~ 3.2G 
  unsigned int start_offset;
  unsigned int end_offset;

  // All exons are marked with the linear location in the entire referenced genome, usually 0 ~ 3.2G
  // This marks the end of the list: exon_ends [total_number_of_exons] = 0
  // It shouldn't be equal to 0, should it be?

  unsigned int exon_starts [MAX_EXONS_PER_GENE];
  unsigned int exon_ends [MAX_EXONS_PER_GENE];
} gene_t;


typedef struct{
	unsigned int memory_block_size;
	unsigned int start_base_offset;
	unsigned int start_point;
	unsigned int length;
	unsigned char * values;
	unsigned int values_bytes;
} gene_value_index_t;



struct gehash_bucket {
	int current_items;
	int space_size;
	union
	{
		short * new_item_keys;
		gehash_key_t * item_keys;
	};
	gehash_data_t * item_values;
};

typedef struct {
	int version_number;
	unsigned long long int current_items;
	int buckets_number;
	char is_small_table;
	struct gehash_bucket * buckets;
	int index_gap;
} gehash_t;


typedef struct {
	gene_vote_number_t max_vote;
	gehash_data_t max_position;
	gene_quality_score_t max_quality;
	gene_vote_number_t max_indel_recorder[MAX_INDEL_TOLERANCE*3];
	gene_vote_number_t * max_tmp_indel_recorder;
	short max_mask;
	gene_vote_number_t noninformative_subreads;

        unsigned short items[GENE_VOTE_TABLE_SIZE];
        unsigned int pos [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
        gene_vote_number_t votes [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
        gene_quality_score_t quality [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short masks [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short last_offset [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	gene_vote_number_t indel_recorder [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE][MAX_INDEL_TOLERANCE*3];
	char current_indel_cursor[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	char toli[GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];

	#ifdef MAKE_FOR_EXON
	short coverage_start [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short coverage_end [GENE_VOTE_TABLE_SIZE][GENE_VOTE_SPACE];
	short max_coverage_start;
	short max_coverage_end;
	//#warning Switch "MAKE_FOR_EXON" is turned on. It may cost more time. Do not turn it on unless you want to detect junction reads.
	#endif
} gene_vote_t ;

typedef struct{
	unsigned int pos;
	char len;
} indel_record_t;

typedef struct{
	int count;
	int support;
} indel_result_t;

typedef struct{
	unsigned char best_len;
	unsigned int offsets [BEXT_RESULT_LIMIT];
	unsigned char is_reverse [BEXT_RESULT_LIMIT];
} gene_best_record_t;



typedef struct{
	int max_len;
	unsigned int * max_positions;
	unsigned char * is_counterpart;
	gene_vote_number_t * max_votes;
	gene_quality_score_t * max_quality;
	gene_quality_score_t * max_final_quality;
	short * masks;
	char * max_indel_recorder;
	char * span_coverage;
#ifdef REPORT_ALL_THE_BEST
	gene_best_record_t * best_records;
#endif
	char max_indel_tolerance;
	short indel_recorder_length;

	unsigned char *repeated_regions;

} gene_allvote_t;


typedef struct{
	int total_offsets;
        char *read_names;
        unsigned int *read_offsets;
	HashTable * read_name_to_index;
} gene_offset_t;


#define EXON_BUFFER_SIZE 3000 

struct thread_input_buffer {
	char read_names [EXON_BUFFER_SIZE][121];
	char read [EXON_BUFFER_SIZE][1201];
	char quality [EXON_BUFFER_SIZE][1201];
	int rl[EXON_BUFFER_SIZE];
	int write_pointer;
	int read_pointer;

	unsigned int read_id[EXON_BUFFER_SIZE];

};


typedef struct {
	char filename [300];
	int space_type ;
	int file_type ;
	FILE * input_fp;
	unsigned int read_chunk_start;
} gene_input_t;


typedef struct{
	unsigned int small_key;
	unsigned int big_key;
} paired_exon_key;

typedef struct{
	unsigned int supporting_reads;
	char is_fusion;
	char big_pos_neg;
	char small_pos_neg;
} fusion_record;

double miltime();

typedef struct{
	char chromosome_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned long known_length;
} chromosome_t;

typedef struct{
	unsigned char record_type;
	unsigned char  mapping_quality;
	unsigned short read_pos;
	unsigned short read_len;
	unsigned short flags;
	unsigned int read_number;
	unsigned int pos;
	char strand;	// 0 = positive, 1 = negative
} base_block_temp_read_t;

typedef struct{
	unsigned char record_type;
	unsigned int pos;
	short type;	
} VCF_temp_read_t;


struct explorer_section_t
{
	unsigned int start_pos;

	short read_pos_start;
	short read_pos_end;

	char is_neg_strand;
	short indels;
	short indel_pos;

	short all_indel_poses[10];
	short all_indels[10];
};

struct explorer_record_t
{
	struct explorer_section_t cigar_record[6];
	short b_search_tail;
};

FILE * f_subr_open(const char * fname, const char * mode);
#define abs(a) 	  ((a)>=0?(a):-(a))
#define max(a,b)  ((a)<(b)?(b):(a))
#define min(a,b)  ((a)>(b)?(b):(a))


#endif
