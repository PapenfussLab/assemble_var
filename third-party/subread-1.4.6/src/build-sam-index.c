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
  
  
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef WINDOWS
#define __i686__
#include <winsock2.h>
#else
#include <arpa/inet.h>
#endif

#include <unistd.h>


#include "subread.h"
#include "sublog.h"
#include "core.h"
#include "input-files.h"
#include "sorted-hashtable.h"

#include "core-indel.h"
#include "core-junction.h"

#define BASE_BLOCK_LENGTH_INDEX 15000000
#define BUCKET_SIZE_INDEX 10000

#define ntohll(x) ( ( (unsigned long long)(ntohl( (unsigned int )(((x) << 32) >> 32) ))<< 32) | ntohl( ((unsigned int)((x) >> 32)) ) )
#define htonll(x) ntohll(x)

char temp_file_prefix[MAX_FILE_NAME_LENGTH];

struct read_position_info
{
	unsigned long long file_offset;
	unsigned int chro_offset;
	unsigned short flags;
	unsigned short section_length;
	char cigar[EXON_MAX_CIGAR_LEN];
};


// BASE_NUMBER_CURVE_POINTS should be (2+4*n) where n is an integer.
#define BASE_NUMBER_CURVE_POINTS 6

struct bucket_info
{
	char chromosome_name[MAX_CHROMOSOME_NAME_LEN];
	unsigned long long int L2index_offset;
	unsigned int chromosome_offset;
	unsigned int bases_in_bucket;
	unsigned int reads_in_bucket;
	unsigned short base_number_curve[BASE_NUMBER_CURVE_POINTS]; 
};

unsigned int calculate_read_span(char * cigar)
{
	int xk0 = 0;
	unsigned int x = 0, ret = 0;
	while (1)
	{
		char nch = cigar[xk0];
		if(!nch)break;
		if(isdigit(nch))
			x = x*10+nch-'0';
		else
		{
			if(nch!='I') ret += x;
			x=0;
		}

		xk0++;
	}

	return ret;
}



unsigned int transform_pillup_to_index(char * chromosome_name, unsigned int block_start, char * temp_file_name, FILE * L1index_fp, FILE * L2index_fp)
{
	FILE * pile_fp = f_subr_open(temp_file_name,"rb");
	if(!pile_fp) return 0;

	int buckets_in_block = BASE_BLOCK_LENGTH_INDEX / BUCKET_SIZE_INDEX;
	struct read_position_info ** buckets = malloc(sizeof(struct read_position_info *) * buckets_in_block);
	unsigned int * bucket_counts = malloc(sizeof(int) * buckets_in_block);
	unsigned int * bucket_sizes = malloc(sizeof(int) * buckets_in_block);
	unsigned int * bucket_bases = malloc(sizeof(int) * buckets_in_block);
	unsigned int * bucket_curve_bases = malloc(sizeof(int) * buckets_in_block * BASE_NUMBER_CURVE_POINTS);
	unsigned int total_reads= 0;

	memset(buckets, 0, sizeof(struct read_position_info *) * buckets_in_block);
	memset(bucket_curve_bases, 0, sizeof(int) * buckets_in_block * BASE_NUMBER_CURVE_POINTS);
	memset(bucket_counts, 0, sizeof(int) * buckets_in_block);
	memset(bucket_sizes, 0, sizeof(int) * buckets_in_block);


	while(!feof(pile_fp))
	{
		struct read_position_info record;
		fread(&record, sizeof(record), 1, pile_fp);
		unsigned int read_chro_offset = ntohl(record.chro_offset), window_pointer; 
		unsigned int read_span = calculate_read_span(record.cigar);

		for(window_pointer = read_chro_offset; window_pointer < read_chro_offset + read_span; window_pointer += BUCKET_SIZE_INDEX )
		{
			if(window_pointer<block_start)continue;
			int bucket_no = (window_pointer-block_start) / BUCKET_SIZE_INDEX;
			if(bucket_no >=buckets_in_block) continue;
			int bucket_curve_offset = (int)((window_pointer-block_start)/(BUCKET_SIZE_INDEX*1./BASE_NUMBER_CURVE_POINTS) - bucket_no * BASE_NUMBER_CURVE_POINTS);

			struct read_position_info * bucket = buckets[bucket_no];
			if(!bucket)
			{
				bucket = malloc(sizeof(struct read_position_info) * 10);
				buckets[bucket_no] = bucket;
				bucket_counts[bucket_no] = 0;
				bucket_sizes[bucket_no] = 10;
			}

			if( bucket_counts[bucket_no] == bucket_sizes[bucket_no])
			{
				bucket = realloc(bucket,  sizeof(struct read_position_info)*bucket_sizes[bucket_no]*1.5);
				buckets[bucket_no] = bucket;
				bucket_sizes[bucket_no] *= 1.5;
			}

			struct read_position_info * item = &(bucket[bucket_counts[bucket_no]]);
			memcpy(item , &record , sizeof(struct read_position_info));
			bucket_counts[bucket_no]++;
			bucket_bases[bucket_no]+= ntohs(record.section_length);
			bucket_curve_bases[bucket_no * BASE_NUMBER_CURVE_POINTS + bucket_curve_offset] += ntohs(record.section_length); 
		}
		total_reads ++;
	}

	int xk1;
	for(xk1 = 0; xk1 < buckets_in_block; xk1++)
	{
		unsigned long long int L2index_offset = ftello(L2index_fp);
		struct bucket_info bucket_record;
		strcpy(bucket_record.chromosome_name, chromosome_name);
		bucket_record.chromosome_offset = htonl(block_start+xk1 * BUCKET_SIZE_INDEX);
		bucket_record.reads_in_bucket = htonl(bucket_counts[xk1]);
		bucket_record.bases_in_bucket = htonl(bucket_bases[xk1]);
		bucket_record.L2index_offset = htonll(L2index_offset);

		int xk2;
		
		for(xk2 = xk1 * BASE_NUMBER_CURVE_POINTS; xk2< xk1 * BASE_NUMBER_CURVE_POINTS +  BASE_NUMBER_CURVE_POINTS; xk2++)
		{
			if(ntohl(bucket_record.bases_in_bucket >0))
			{
				bucket_record.base_number_curve[xk2 - xk1 * BASE_NUMBER_CURVE_POINTS] = htons((unsigned short)(bucket_curve_bases[xk2] * 30000./ ntohl(bucket_record.bases_in_bucket)));
				//printf("RV=%d; CV=%d; TOTAL=%d\n", bucket_curve_bases[xk2]  , ntohs(bucket_record.base_number_curve[xk2 - xk1 * BASE_NUMBER_CURVE_POINTS]), ntohl( bucket_record.bases_in_bucket));
			}
			else	bucket_record.base_number_curve[xk2 - xk1 * BASE_NUMBER_CURVE_POINTS] = 0;
		}

		fwrite(&bucket_record, sizeof(struct bucket_info), 1, L1index_fp);
		fwrite(buckets[xk1] , sizeof(struct read_position_info) * bucket_counts[xk1],1, L2index_fp);
	}

	fclose(pile_fp);
	unlink(temp_file_name);

	free(buckets);
	free(bucket_curve_bases);
	free(bucket_counts);
	free(bucket_sizes);

	return total_reads;
}

int finalise_sam_index(HashTable * chromosome_size_table, char * output_file_prefix)
{
	KeyValuePair * cursor;
	char temp_file_name[MAX_FILE_NAME_LENGTH];
	sprintf(temp_file_name, "%s.L1i",output_file_prefix);
	FILE * L1index_fp = f_subr_open(temp_file_name,"wb");
	sprintf(temp_file_name, "%s.L2i",output_file_prefix);
	FILE * L2index_fp = f_subr_open(temp_file_name,"wb");
	int bucket;
	for(bucket=0; bucket<chromosome_size_table -> numOfBuckets; bucket++)
	{
		cursor = chromosome_size_table -> bucketArray[bucket];

		while(1)
		{
			if(!cursor) break;
			char * chro_name = (char *)cursor -> key;
			unsigned int chro_max_length =cursor ->value - NULL - 1;

			unsigned int chro_offset = 0;


			SUBREADprintf(" === %s : 0 ~ %u === \n", chro_name , chro_max_length);

			for(chro_offset = 0; chro_offset < chro_max_length; chro_offset += BASE_BLOCK_LENGTH_INDEX)
			{

				sprintf(temp_file_name,"%s@%s-%04u.bin", temp_file_prefix, chro_name , chro_offset / BASE_BLOCK_LENGTH_INDEX );

				unsigned int total_reads = transform_pillup_to_index(chro_name, chro_offset, temp_file_name, L1index_fp, L2index_fp);
				SUBREADprintf("%s has %d reads\n", temp_file_name, total_reads);
			}
			free(chro_name);
			cursor = cursor->next;
		}
	}

	fclose(L1index_fp);
	fclose(L2index_fp);
	return 0;
}

int write_read_pos(char * chro_name, unsigned int chro_offset_anchor, unsigned short section_chro_length, short flags, unsigned long long file_offset, char * cigar_str, unsigned int read_span,  HashTable *pileup_fp_table)
{
	char temp_file_name[MAX_FILE_NAME_LENGTH];
	unsigned int window_point;
	unsigned int last_window_no = 0xffffffff;

	for(window_point = chro_offset_anchor; window_point < chro_offset_anchor+read_span; window_point += BUCKET_SIZE_INDEX)
	{
		unsigned int window_no = window_point / BASE_BLOCK_LENGTH_INDEX;

		if(window_no != last_window_no)
		{

			sprintf(temp_file_name,"%s@%s-%04u.bin", temp_file_prefix, chro_name , window_point / BASE_BLOCK_LENGTH_INDEX );

			FILE * pileup_fp = get_temp_file_pointer(temp_file_name, pileup_fp_table); 
			struct read_position_info record;

			record.chro_offset = htonl(chro_offset_anchor);
			record.file_offset = htonll(file_offset);
			record.section_length = htons(section_chro_length);
			record.flags = htons(flags);
			strcpy(record.cigar, cigar_str);

			fwrite(&record, sizeof(record), 1, pileup_fp);

			last_window_no = window_no;
		}
	}
	return 0;
}	


// returns 0 if there is no more sections; or 1.
int next_read_section(char * cigar, unsigned int * chro_offset, unsigned short * read_offset, unsigned short * section_chro_length, int * cigar_cursor)
{
	int x=0, local_cigar_cursor = 0;
	int nch;
	unsigned int ret=0;
	unsigned int all_chro_offset = 0;
	unsigned short all_read_offset = 0;
	unsigned int all_section_chro_length = 0;
	

	while(1)
	{
		nch = cigar[local_cigar_cursor];
		if(isdigit(nch))
		{
			x = x*10+(nch-'0');
		}
		else{
			if(nch=='M' || nch == 'S')
			{
				(*chro_offset) = all_chro_offset;
				(*read_offset) = all_read_offset;
				(*section_chro_length) = all_section_chro_length;
			}

			if(nch=='M'||nch=='I' || nch=='S')
				all_read_offset += x;
			if(nch=='D'||nch=='N' || nch=='M' || nch=='S')
				all_chro_offset += x;
			if(nch=='M'||nch=='S'||nch=='D')
				all_section_chro_length += x;
			if(nch=='N') all_section_chro_length = 0;

			x=0;

			if((local_cigar_cursor> (*cigar_cursor) && (nch=='N' || nch == 'I' || nch == 'D')) || cigar[local_cigar_cursor+1] == 0)
			{
				(*cigar_cursor) = local_cigar_cursor+1;
				return !(cigar[local_cigar_cursor+1] == 0);
			}
		
		}

		local_cigar_cursor++;
		if((*cigar_cursor) == local_cigar_cursor)all_section_chro_length = 0;

		if( cigar[local_cigar_cursor]==0) return 0;

	}

	return ret;
}
int transfer_SAM_to_position_table(char * sam_file)
{
	char linebuf[max(2*MAX_READ_LENGTH+300,3000)];
	int allreads=0,mapped=0;
	char mate_chro[MAX_CHROMOSOME_NAME_LEN+1];
	unsigned int mate_pos;
	int flags_mate;
	char cigar_mate[EXON_MAX_CIGAR_LEN+1];
	gene_input_t input_file;

	HashTable * local_reassembly_pileup_files;
	HashTable * chromosome_size_table;


	chromosome_size_table = HashTableCreate(100);
	HashTableSetDeallocationFunctions(chromosome_size_table, NULL, NULL);
	HashTableSetKeyComparisonFunction(chromosome_size_table, my_strcmp);
	HashTableSetHashFunction(chromosome_size_table, HashTableStringHashFunction);

	local_reassembly_pileup_files = HashTableCreate(100);
	HashTableSetDeallocationFunctions(local_reassembly_pileup_files, NULL, NULL);
	HashTableSetKeyComparisonFunction(local_reassembly_pileup_files, my_strcmp);
	HashTableSetHashFunction(local_reassembly_pileup_files ,HashTableStringHashFunction);

	geinput_open_sam(sam_file, &input_file,0);

	while(1)
	{
		unsigned int read_pos =0xffffffff;
		char read_name[MAX_READ_NAME_LEN+1];
		int flags;
		char chro_name[MAX_CHROMOSOME_NAME_LEN+1];
		unsigned int pos, pair_dist;
		char cigar[EXON_MAX_CIGAR_LEN+1];
		char read_text[MAX_READ_LENGTH+1];
		char qual_text[MAX_READ_LENGTH+1];
		int read_len, is_repeated, mapping_quality, is_anchor_certain=1;
		int pos_delta = 0;
		int is_paired_end_reads = 0;

		unsigned long long int file_offset = ftello(input_file.input_fp);
		if(feof(input_file.input_fp))break;

		read_text[0]=0;
		qual_text[0]=0;
		geinput_readline(&input_file, linebuf,0);
		int res = parse_SAM_line(linebuf, read_name,& flags, chro_name, & pos, cigar, & mapping_quality, & pair_dist, read_text , qual_text, & read_len, & is_repeated);
		int cigar_cursor = 0;
		int firstM = 1,xx=0;
		is_paired_end_reads = flags &1;
		if(res==0)
		{
			for(; cigar[cigar_cursor]; cigar_cursor++)
			{
				char nch = cigar[cigar_cursor]; 
				if(nch>='0'&&nch<='9')
				{
					xx=xx*10+(nch-'0');
				}else
				{
					if(nch=='M') firstM=0;
					else if(nch=='S')
					{
						if(firstM)
							pos_delta = xx;
					}
					xx=0;
				}
			}

			pos -= pos_delta;
		}

		if(res == 1) {chro_name[0]='*'; chro_name[1]=0;}
		//printf("MAPPED=%d\n", res);

		if(res == 0)	// mapped
		{

			read_pos = pos - 1;
			mapped++;
		}
		else if(res == 1 && is_paired_end_reads)	// unmapped
		{
			is_anchor_certain=0;
			if(mate_chro[0])
			{
				if(mate_chro[0]!='*')
				{
					read_pos = mate_pos ;
					strcpy(chro_name, mate_chro);
				}
				//printf("RECOVERED 1: %u - %s ; LEN=%d ; YOU_ARE_FIRST=%d\n%s\n%s\n", read_anchor_position,  read_name, read_len, flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR, read_text, qual_text);
			}
			else
			{
				char read_text_null[MAX_READ_LENGTH+1];
				char qual_text_null[MAX_READ_LENGTH+1];

				geinput_readline_back(&input_file, linebuf);
				res = parse_SAM_line(linebuf, read_name,& flags_mate, mate_chro, & mate_pos, cigar_mate, & mapping_quality, & pair_dist, read_text_null , qual_text_null, & read_len, & is_repeated);
				if(res==0)
				{
					read_pos = mate_pos - 1;
					strcpy(chro_name, mate_chro);
				}
				//printf("RECOVERED 2: %u - %s ; LEN=%d ; YOU_ARE_FIRST=%d\n%s\n%s\n", read_anchor_position,  read_name, read_len, flags_mate & SAM_FLAG_FIRST_READ_IN_PAIR, read_text, qual_text);
			}
		}

		if(read_pos<0xffff0000)
		{

			unsigned int read_span = calculate_read_span(cigar);

			write_read_pos(chro_name, read_pos, read_len, flags, file_offset, cigar, read_span, local_reassembly_pileup_files);
			unsigned int old_max_pos = HashTableGet(chromosome_size_table, chro_name) - NULL;
			if(old_max_pos==0)
			{
				char * chro_name_new = malloc(strlen(chro_name)+1);
				strcpy(chro_name_new , chro_name);
				
				HashTablePut(chromosome_size_table, chro_name_new, NULL+read_pos + read_span+1);
			}
			else if(read_pos + read_span +1 > old_max_pos )
				HashTablePutReplace(chromosome_size_table, chro_name, NULL+read_pos + read_span+1, 0);
		}

		if(allreads %2==1)
			mate_chro[0] = 0;
		else
		{
			strcpy(mate_chro, chro_name);
			mate_pos = pos - 1;
			flags_mate = flags;
			
		}
		allreads++;

	}
	SUBREADprintf("Processed %d reads; %d mapped.\n", allreads, mapped);

	destroy_pileup_table(local_reassembly_pileup_files);


	finalise_sam_index(chromosome_size_table, sam_file);

	HashTableDestroy(chromosome_size_table);

	return 0;
}


#ifdef MAKE_STANDALONE
int main(int argc , char ** argv)
{
	//progress_report_callback = NULL;
#else
int samindex_main(int argc , char ** argv)
{
#endif

	char c;
	char in_SAM_file[MAX_FILE_NAME_LENGTH+1];
	optind = 1;
	opterr = 1;
	optopt = 63;

	while ((c = getopt(argc, argv, "i:?")) != -1)
	{
		switch(c)
		{
			case 'i':
				strncpy(in_SAM_file, optarg, MAX_FILE_NAME_LENGTH);
			break;

		}
	}
	sprintf(temp_file_prefix, "./index-temp-sum-%06u-%06u", getpid(),rand());

	int ret = 0;
	ret = transfer_SAM_to_position_table(in_SAM_file);
	return ret;
}
