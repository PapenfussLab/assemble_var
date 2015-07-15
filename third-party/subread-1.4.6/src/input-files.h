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
  
  
#ifndef __INPUT_FILES_H_
#define __INPUT_FILES_H_

#include "subread.h"
#include "hashtable.h"
#include "core-indel.h"

#define GENE_SPACE_BASE 1
#define GENE_SPACE_COLOR 2

#define GENE_INPUT_PLAIN 0
#define GENE_INPUT_FASTQ 1
#define GENE_INPUT_FASTA 2

#define GENE_INPUT_SAM_SINGLE   93
#define GENE_INPUT_SAM_PAIR_1   94
#define GENE_INPUT_SAM_PAIR_2   95



#define FILE_TYPE_SAM     50
#define FILE_TYPE_BAM     500
#define FILE_TYPE_FAST_   100
#define FILE_TYPE_FASTQ   105
#define FILE_TYPE_FASTA   110
#define FILE_TYPE_GZIP_FAST_   1000
#define FILE_TYPE_GZIP_FASTQ   1105
#define FILE_TYPE_GZIP_FASTA   1110
#define FILE_TYPE_UNKNOWN 999
#define FILE_TYPE_EMPTY   999990
#define FILE_TYPE_NONEXIST 999999



#include <stdlib.h>
#include <stdio.h>



#define SAM_SORT_BLOCKS 229
#define SAM_SORT_BLOCK_SIZE 512333303LLU
//#define SAM_SORT_BLOCK_SIZE 11123333LLU

typedef struct
{
	unsigned long long int output_file_size;
	unsigned long long int current_chunk_size;
	unsigned int current_chunk;
	unsigned long long int written_reads;
	unsigned long long int unpaired_reads;
	FILE * current_block_fp_array [SAM_SORT_BLOCKS];
	FILE * all_chunks_header_fp;

	FILE * out_fp;
	char tmp_path[MAX_FILE_NAME_LENGTH];
} SAM_sort_writer;


void fastq_64_to_33(char * qs);

int chars2color(char c1, char c2);

int genekey2color(char last_base,char * key);

// Convert a string key into an integer key
int genekey2int(char key [], int space_type);

// Open a read file. This function automatically decides its type.
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open(const char * filename, gene_input_t * input);

// Open a sam file. Parameter half_no indicates which read of the pair is concerned.
// half_no = 1: the first read; half_no = 2: the last read; half_no = 0: single-end
// Return 0 if successfully opened, or 1 if error occurred
int geinput_open_sam(const char * filename, gene_input_t * input, int half_no);

// Read a line from the input and reward the pointer to the last position.
int geinput_readline_back(gene_input_t * input, char * linebuffer) ;

// Get the next read from the input file
// Return the length of this read or -1 if EOF. 
// The memory space for read_string must be at least 512 bytes.
int geinput_next_read(gene_input_t * input, char * read_name, char * read_string, char * quality_string);
int geinput_next_read_sam(gene_input_t * input, char * read_name, char * read_string, char * quality_string, gene_offset_t* offsets, unsigned int * pos, int * mapping_quality, int * mapping_flags, int need_reversed);
int geinput_next_read_trim(gene_input_t * input, char * read_name, char * read_string, char * quality_string, short trim_5, short trim_3, int * is_secondary);

void geinput_jump_read(gene_input_t * input);

// Close the input file
void geinput_close(gene_input_t * input);

// return the next ATGC char from a input file.
// return 0 if this read segment reaches the end; return -1 if EOF; return -2 if error
int geinput_next_char(gene_input_t * input);

// line buffer has to be at least 300 bytes
// it returns the length of reading
int geinput_readline(gene_input_t * input, char * linebuffer, int conv_to_upper) ;



// read a line into the buff,
// the line should not be longer than 300 chars or the remaining part will be discarded.
// therefore the buff has to be at least 300 chars.
int read_line(int max_len, FILE * fp, char * buff, int conv_to_upper);

// count the number of reads in a flie
double guess_reads_density(char * fname, int is_sam) ;

// guess the size of the chromosome lib
// return the number of bases, or (-index-1) if the file at the index is not found.
long long int guess_gene_bases(char ** files, int file_number);


void reverse_read(char * ReadString, int Length, int space_type);

void reverse_quality(char * QualtyString, int Length);

unsigned int read_numbers(gene_input_t * input);

//This function returns 0 if the line is a mapped read; -1 if the line is in a wrong format and 1 if the read is unmapped.
int parse_SAM_line(char * sam_line, char * read_name, int * flags, char * chro, unsigned int * pos, char * cigar, int * mapping_quality, unsigned int * pair_dist, char * sequence , char * quality_string, int * rl, int * repeated);

#define reverse_char(c)	((c)=='A'?'T':((c)=='G'?'C':((c)=='C'?'G':'A')))

int find_subread_end(int len, int  TOTAL_SUBREADS,int subread) ;

int break_SAM_file(char * in_SAM_file, int is_BAM, char * temp_file_prefix, unsigned int * real_read_count, int * block_no, chromosome_t * known_chromosomes, int is_sequence_needed, int base_ignored_head_tail, gene_value_index_t *array_index, gene_offset_t * offsets, unsigned long long int * all_Mapped_bases , HashTable * event_table_ptr, char * VCF_file);

int get_known_chromosomes(char * in_SAM_file, chromosome_t * known_chromosomes);


int load_exon_annotation(char * annotation_file_name, gene_t ** output_genes, gene_offset_t* offsets );

int is_in_exon_annotations(gene_t *output_genes, unsigned int offset, int is_start);

int does_file_exist (char * filename);

double guess_reads_density_format(char * fname, int is_sam, int * min_phred, int * max_phred);

FILE * get_temp_file_pointer(char *temp_file_name, HashTable* fp_table);

void write_read_block_file(FILE *temp_fp , unsigned int read_number, char *read_name, int flags, char * chro, unsigned int pos, char *cigar, int mapping_quality, char *sequence , char *quality_string, int rl , int is_sequence_needed, char strand, unsigned short read_pos, unsigned short read_len);

int get_read_block(char *chro, unsigned int pos, char *temp_file_suffix, chromosome_t *known_chromosomes, unsigned int * max_base_position);
int my_strcmp(const void * s1, const void * s2);

void destroy_cigar_event_table(HashTable * event_table);


int is_SAM_unsorted(char * SAM_line, char * tmp_read_name, short * tmp_flag, unsigned long long int read_no);
int sort_SAM_add_line(SAM_sort_writer * writer, char * SAM_line, int line_len);
void sort_SAM_finalise(SAM_sort_writer * writer);
int sort_SAM_create(SAM_sort_writer * writer, char * output_file, char * tmp_path);
void colorread2base(char * read_buffer, int read_len);

int warning_file_type(char * fname, int expected_type);
char color2char(char clr, char c1);

int is_certainly_bam_file(char * fname, int * is_firstread_PE);

unsigned long long int sort_SAM_hash(char * str);

char * fgets_noempty(char * buf, int maxlen, FILE * fp);

char * gzgets_noempty(void * fp, char * buf, int maxlen);
#endif
