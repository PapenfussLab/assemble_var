#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "subread.h" 
#include "core.h" 
#include "HelperFunctions.h"
#include "sambam-file.h" 
#include "gene-algorithms.h" 
#include "input-files.h" 
#include "hashtable.h"

#define COVERAGE_MAX_INT 254
unsigned long long all_counted;
typedef unsigned int coverage_bin_entry_t;
int is_BAM_input = 0;
char input_file_name[300];
char output_file_name[300];
HashTable * cov_bin_table;


static struct option cov_calc_long_options[] =
{
	{"primary",no_argument, 0, 0},
	{0, 0, 0, 0}
};


void calcCount_usage()
{
	SUBREADprintf("\ncoveageCount v%s\n\n", SUBREAD_VERSION);
	SUBREADprintf("Counting the coverage of mapped reads at each location on the entire reference genome.\n\n");
	SUBREADprintf("./ncoveageCount -i <sam_bam_input> -o <output_prefix>\n\n");
}

void add_chro(char *sam_h)
{
	char *chro_name = malloc(200);
	unsigned int chro_len = 0;

	char nch;
	int cur = 0, tabs=0, state = 0, txtcur = 0;

	chro_name[0]=0;
	while(1){
		nch = sam_h[cur++];
		if(!nch || nch == '\n')break;

		if(nch == '\t')
		{
			txtcur = 0;
			tabs ++;
			state = 0;
		}
		else if(nch == ':')
		{
			state ++;
		}
		else
		{
			if(state == 1 && tabs == 1)
			{
				chro_name[txtcur++]=nch;
				chro_name[txtcur]=0;
			}
			else if(state == 1 && tabs ==2)
				chro_len = chro_len*10 + (nch-'0');
		}
	}

	if(chro_name[0]==0 || chro_len<1)
	{
		SUBREADprintf("ERROR: incorrect SAM format: %s\n", sam_h);
	}

	void ** bin_entry = malloc(sizeof(void *)*2);
	coverage_bin_entry_t * chro_bin = calloc(sizeof(coverage_bin_entry_t) , chro_len);
	if(!chro_bin)
	{
		SUBREADprintf("ERROR: cannot allocate the memory block. You need at least 4GB of memory,\n");
	}
	bin_entry[0] = (void *)chro_bin;
	bin_entry[1] = (void *)(NULL + chro_len);
	HashTablePut(cov_bin_table, chro_name, bin_entry);
	SUBREADprintf("Added a new chromosome : %s [%u]\n", chro_name, chro_len);
}

void get_read_info(char * fl, char * chro, unsigned int * pos , char * cigar, int *flags){
	char * tmp_tok = NULL;

	char * tmp_res = strtok_r(fl, "\t", &tmp_tok);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok);
	(*flags) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // chro
	strcpy(chro, tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // pos
	(*pos) = atoi(tmp_res);

	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // qual
	tmp_res = strtok_r(NULL, "\t", &tmp_tok); // cigar
	strcpy(cigar, tmp_res);
}

int covCalc()
{

	cov_bin_table = HashTableCreate(200);
	HashTableSetHashFunction(cov_bin_table, fc_chro_hash);
	HashTableSetKeyComparisonFunction(cov_bin_table , fc_strcmp_chro);
	HashTableSetDeallocationFunctions(cov_bin_table , free, free);


	SamBam_FILE * in_fp = SamBam_fopen(input_file_name, is_BAM_input?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	char * line_buffer = malloc(3000);
	
	while(1)
	{
		char * is_ret = SamBam_fgets(in_fp, line_buffer, 2999, 1);
		if(!is_ret) break;
		if(line_buffer[0]=='@'){
			if(strstr(line_buffer,"@SQ\t"))
				add_chro(line_buffer);
		}
		else
		{
			unsigned int Staring_Points[6];
			unsigned short Section_Lengths[6];
	
			int flags=0, x1;
			char cigar_str[200];
			char chro[200];
			unsigned int pos = 0;
			cigar_str[0]=0;
			chro[0]=0;

			get_read_info(line_buffer, chro, &pos, cigar_str, &flags);

			if(flags & 4) continue;

			void ** bin_entry = HashTableGet(cov_bin_table, chro);
			if(NULL == bin_entry)
			{
				SUBREADprintf("ERROR: The chromosome name is not in header:%s\n", chro);
			}

			coverage_bin_entry_t * chrbin = (coverage_bin_entry_t*) bin_entry[0];
			unsigned int chrlen = (void *)( bin_entry[0]) - NULL;
			int cigar_sections = RSubread_parse_CIGAR_string(cigar_str, Staring_Points, Section_Lengths);
			for(x1 = 0; x1 < cigar_sections; x1++)
			{
				unsigned int x2;
				//printf("%s: %u-%u\n", chro, Staring_Points[x1], Staring_Points[x1]+Section_Lengths[x1]);
				for(x2 = Staring_Points[x1]; x2<Staring_Points[x1]+Section_Lengths[x1]; x2++)
				{
					if(pos+x2 < chrlen)
					{
						if(chrbin[x2 + pos] <= COVERAGE_MAX_INT)chrbin[x2 + pos] ++;
						all_counted ++;
						if(all_counted % 10000000 == 0)
						{
							SUBREADprintf("Processed %llu bases.\n", all_counted);
						}
					}
				}
			}
		}
	}

	free(line_buffer);

	SamBam_fclose(in_fp);


	SUBREADprintf("Processed totally %llu bases.\nNow write results.\n", all_counted);

	int bucket;
	KeyValuePair *cursor;
	for(bucket=0; bucket < cov_bin_table  -> numOfBuckets; bucket++)
	{
		cursor = cov_bin_table -> bucketArray[bucket];
		while (1)
		{
			if (!cursor) break;
			coverage_bin_entry_t * chrbin = (coverage_bin_entry_t*)(((void **) cursor -> value)[0]);
			unsigned int chrlen = (((void **) cursor -> value)[1]) - NULL; 
			char * chro = (char *)(cursor -> key);
			char out_name[340];
			sprintf(out_name,"%s-%s.bin", output_file_name, chro);

			FILE * fpo = fopen(out_name,"w");
			fwrite(chrbin, sizeof(coverage_bin_entry_t), chrlen, fpo);
			fclose(fpo);
			free(chrbin);

			SUBREADprintf("Wrote bin for %s\n", chro);
			cursor = cursor->next;
		}	
	}

	HashTableDestroy(cov_bin_table);

	SUBREADprintf("Calculation finished.\n");
	return 0;
}



#ifdef MAKE_STANDALONE
int main(int argc, char ** argv)
#else
int cov_calc_main(int argc, char ** argv)
#endif
{
	int ret = 0;
	int c=0;
	int option_index=0;
	input_file_name[0]=0;
	output_file_name[0]=0;
	is_BAM_input=0;
	all_counted = 0;

	optind=0;
	opterr=1;
	optopt=63;

	while ((c = getopt_long (argc, argv, "Bi:o:?", cov_calc_long_options, &option_index)) != -1)
		switch(c)
		{
			case 'i':
				strcpy(input_file_name, optarg);
			break;
			case 'o':
				strcpy(output_file_name, optarg);
			break;

			case '?':
			default :
				calcCount_usage();
				return -1;
	
		}
	

	if((!output_file_name[0])||(!input_file_name[0]))
	{
		calcCount_usage();
		return 0;
	}

	int is_bam = is_certainly_bam_file(input_file_name, NULL);

	if(1==is_bam) is_BAM_input = 1;
	else if(is_bam < 0)
	{
		ret = -1;
		SUBREADprintf("Unable to open input file '%s' or the input file is empty!\n", input_file_name);
	}

	ret = ret || covCalc();

	return ret;
}
