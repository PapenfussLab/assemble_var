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

/***************************************************************

   The ASCII Art used in this file was generated using FIGlet and
   the big font, contributed by Glenn Chappell to FIGlet.
  
  ***************************************************************/
  
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/stat.h>
#include <ctype.h>



#include "subread.h"
#include "sublog.h"
#include "core.h"
#include "input-files.h"
#include "sorted-hashtable.h"

#include "core-indel.h"
#include "core-junction.h"

static struct option long_options[] =
{
	{"memory-optimisation",  required_argument, 0, 0},
	{0, 0, 0, 0}
};

int (*progress_report_callback)(int, int, int);

int is_result_in_PE(alignment_result_t *al)
{
	if(al->Score_H & 0x8000000000000000llu)return 1;
	return 0;
}

void core_version_number(char * program)
{
	SUBREADprintf("\n%s v%s\n\n" , program, SUBREAD_VERSION);
}

void warning_file_limit()
{
	struct rlimit limit_st;
	getrlimit(RLIMIT_NOFILE, & limit_st);

	{
		if(min(limit_st.rlim_cur , limit_st.rlim_max) < 400)
		{
			print_in_box(80,0,0,"WARNING This operation needs to open many files at same time,");
			print_in_box(80,0,0,"	but your OS only allows to open %d files.", min(limit_st.rlim_cur , limit_st.rlim_max));
			print_in_box(80,0,0,"	You can use command 'ulimit -n 500' to raise this limit");
			print_in_box(80,0,0,"	to 500, or the program may crash or terminate unexpectedly.");
			print_in_box(80,0,0,"");
		}
	}
}

void print_in_box(int line_width, int is_boundary, int is_center, char * pattern,...)
{
	va_list args;
	va_start(args , pattern);
	char is_R_linebreak=0, * content, *out_line_buff;

	content= malloc(1000);
	out_line_buff= malloc(1000);
	out_line_buff[0]=0;
	vsprintf(content, pattern, args);
	int is_R_code,x1,content_len = strlen(content), state, txt_len, is_cut = 0, real_lenwidth;

	is_R_code = 1;
	#ifdef MAKE_STANDALONE
		is_R_code = 0;
	#endif

	if(content_len>0&&content[content_len-1]=='\r'){
		content_len--;
		content[content_len] = 0;
		is_R_linebreak = 1;
	}

	if(content_len>0&&content[content_len-1]=='\n'){
		content_len--;
		content[content_len] = 0;
	}

	state = 0;
	txt_len = 0;
	real_lenwidth = line_width;
	for(x1 = 0; content [x1]; x1++)
	{
		char nch = content [x1];
		if(nch == CHAR_ESC)
			state = 1;
		if(state){
			real_lenwidth --;
		}else{
			txt_len++;
			
			if(txt_len == 80 - 6)
			{
				is_cut = 1;
			} 
		}

		if(nch == 'm' && state)
			state = 0;
	}

	if(is_cut)
	{
		state = 0;
		txt_len = 0;
		for(x1 = 0; content [x1]; x1++)
		{
			char nch = content [x1];
			if(nch == CHAR_ESC)
				state = 1;
			if(!state){
				txt_len++;
				if(txt_len == 80 - 9)
				{
					strcpy(content+x1, "\x1b[0m ...");
					content_len = line_width - 4;
					content_len = 80 - 4;
					line_width = 80;
					break;
				} 
			}
			if(nch == 'm' && state)
				state = 0;
		}
	}

	if(content_len==0 && is_boundary)
	{
		strcat(out_line_buff,is_boundary==1?"//":"\\\\");
		for(x1=0;x1<line_width-4;x1++)
			strcat(out_line_buff,"=");
		strcat(out_line_buff,is_boundary==1?"\\\\":"//");
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);

		free(content);
		free(out_line_buff);
		return;
	}
	else if(is_boundary)
	{
		int left_stars = (line_width - content_len)/2 - 1;
		int right_stars = line_width - content_len - 2 - left_stars;
		strcat(out_line_buff,is_boundary==1?"//":"\\\\");
		for(x1=0;x1<left_stars-2;x1++) strcat(out_line_buff,"=");
		sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
		sprintf(out_line_buff+strlen(out_line_buff)," %s ", content);
		sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
		for(x1=0;x1<right_stars-2;x1++) strcat(out_line_buff,"=");
		strcat(out_line_buff,is_boundary==1?"\\\\":"//");
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "%s", out_line_buff);

		free(content);
		free(out_line_buff);
		return;
	}

	int right_spaces, left_spaces;	
	if(is_center)
		left_spaces = (line_width - content_len)/2-2;
	else
		left_spaces = 1;

	right_spaces = line_width - 4 - content_len- left_spaces; 

	char spaces[81];
	memset(spaces , ' ', 80);
	spaces[0]='|';
	spaces[1]='|';
	spaces[80]=0;

	//sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"||");
	
	//for(x1=0;x1<left_spaces;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");

	spaces[left_spaces+2] = 0;
	strcat(out_line_buff,spaces);

	if(is_R_code)
	{
		strcat(out_line_buff,content);
	}
	else
	{
		int col1w=-1;
		for(x1=0; content[x1]; x1++)
		{
			if(content[x1]==':')
			{
				col1w=x1;
				break;
			}
		}
		if(col1w>0 && col1w < content_len-1)
		{
			content[col1w+1]=0;
			strcat(out_line_buff,content);
			strcat(out_line_buff," ");
			sprintf(out_line_buff+strlen(out_line_buff),"%c[36m", CHAR_ESC);
			strcat(out_line_buff,content+col1w+2);
			sprintf(out_line_buff+strlen(out_line_buff),"%c[0m", CHAR_ESC);
		}
		else
			strcat(out_line_buff,content);
	}
//	for(x1=0;x1<right_spaces - 1;x1++) sublog_fwrite(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO," ");
	
	memset(spaces , ' ', 80);
	spaces[79]='|';
	spaces[78]='|';
	
	right_spaces = max(1,right_spaces);
	//if(is_R_linebreak)
	//	sprintf(out_line_buff+strlen(out_line_buff)," %c[0m%s%c", CHAR_ESC, spaces + (78 - right_spaces + 1) ,CORE_SOFT_BR_CHAR);
	//else
		sprintf(out_line_buff+strlen(out_line_buff)," %c[0m%s", CHAR_ESC , spaces + (78 - right_spaces + 1));
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, out_line_buff);
	free(out_line_buff);
	free(content);
}




int show_summary(global_context_t * global_context)
{

        if(progress_report_callback)
        {
                long long int all_reads_K = global_context -> all_processed_reads / 1000;
                float mapped_reads_percentage = global_context -> all_mapped_reads * 1./global_context -> all_processed_reads;
                if(global_context->input_reads.is_paired_end_reads) mapped_reads_percentage/=2;
                progress_report_callback(10, 900000, (int) (miltime()-global_context->start_time));
                progress_report_callback(10, 900010, (int) all_reads_K);
                progress_report_callback(10, 900011, (int) (10000.*mapped_reads_percentage));
        }

        print_in_box(80,0,1,"");
        print_in_box(89,0,1,"%c[36mCompleted successfully.%c[0m", CHAR_ESC, CHAR_ESC);
        print_in_box(80,0,1,"");
        print_in_box(80,2,1,"");
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
        print_in_box(80, 1,1,"Summary");
        print_in_box(80, 0,1,"");
        print_in_box(80, 0,0,"         Processed : %llu %s" , global_context -> all_processed_reads, global_context->input_reads.is_paired_end_reads?"fragments":"reads");
        print_in_box(81, 0,0,"            Mapped : %llu %s (%.1f%%%%)", global_context -> all_mapped_reads, global_context->input_reads.is_paired_end_reads?"fragments":"reads" ,  global_context -> all_mapped_reads*100.0 / global_context -> all_processed_reads);
        if(global_context->input_reads.is_paired_end_reads)
                print_in_box(80, 0,0,"  Correctly paired : %llu fragments", global_context -> all_correct_PE_reads);

        if(global_context->config.output_prefix[0])
        {
                if(global_context->config.entry_program_name == CORE_PROGRAM_SUBJUNC)
                        print_in_box(80, 0,0,"         Junctions : %u", global_context -> all_junctions);
                if(global_context->config.do_fusion_detection)
                        print_in_box(80, 0,0,"           Fusions : %u", global_context -> all_fusions);
                print_in_box(80, 0,0,"            Indels : %u", global_context -> all_indels);
        }
        

	if(global_context -> is_phred_warning)
	{
		print_in_box(80, 0,1,"");
		print_in_box(80,0,0, "           WARNING : Phred offset (%d) incorrect?", global_context->config.phred_score_format == FASTQ_PHRED33?33:64);
	}
        print_in_box(80, 0,1,"");
        print_in_box(80, 0,0,"      Running time : %.1f minutes", (miltime()-global_context->start_time)*1./60);
        print_in_box(80, 0,1,"");
        print_in_box(80, 2,1,"http://subread.sourceforge.net/");
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");

        return 0;
}



void show_progress(global_context_t * global_context, thread_context_t * thread_context, unsigned int current_read_no, int task)
{
	if(thread_context&&thread_context->thread_id)
	{
		SUBREADputs("show_progress can only be called by thread#0\n");
		return;
	}

	gene_input_t * ginp1 = thread_context?(thread_context->ginp1):(&global_context->input_reads.first_read_file);
	unsigned long long ginp1_file_pos = ftello(ginp1->input_fp);

	if(task == STEP_VOTING)
	{
		unsigned long long real_read_number = global_context -> all_processed_reads + current_read_no;// * global_context -> config.all_threads;
		if(real_read_number>1000)
			global_context -> input_reads . avg_read_length = (ginp1_file_pos - ginp1->read_chunk_start) * 1./real_read_number ;
	}

	unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
	unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;
	//printf("FS=%llu; AVG=%f; GAR=%llu; CURRENT_NO=%u\n", total_file_size, global_context -> input_reads . avg_read_length , guessed_all_reads, current_read_no);

	unsigned long long current_block_start_file_offset = global_context -> current_circle_start_position_file1;

	unsigned long long guessed_this_chunk_all_reads = (total_file_size - current_block_start_file_offset) / global_context -> input_reads . avg_read_length ;
	if(guessed_this_chunk_all_reads > global_context ->config.reads_per_chunk) guessed_this_chunk_all_reads = global_context ->config.reads_per_chunk;

	unsigned long long guessed_all_reads_before_this_chunk = current_block_start_file_offset / global_context -> input_reads . avg_read_length ;
	unsigned long long reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;//* global_context -> config.all_threads;
	if(task != STEP_VOTING)
	   reads_finished_in_this_chunk = (ginp1_file_pos - current_block_start_file_offset) / global_context -> input_reads . avg_read_length;
	

	unsigned long long finished_steps = guessed_all_reads_before_this_chunk * (global_context -> index_block_number * 6 + 4);
	
	if(task == STEP_VOTING)
		finished_steps += guessed_this_chunk_all_reads * global_context -> current_index_block_number * 6 ;//* global_context -> config.all_threads; 
	if(task >= STEP_ITERATION_ONE)
		finished_steps += guessed_this_chunk_all_reads * global_context -> index_block_number * 6 ;//*global_context -> config.all_threads;
	if(task > STEP_ITERATION_ONE)
		finished_steps += guessed_this_chunk_all_reads ;//* global_context -> config.all_threads;
	if(task > STEP_ITERATION_TWO)
		finished_steps += guessed_this_chunk_all_reads;

	if(task == STEP_VOTING)	finished_steps += reads_finished_in_this_chunk*5*global_context -> config.all_threads;
	if(task > STEP_ITERATION_TWO)
		finished_steps += reads_finished_in_this_chunk * 2;
	else
		finished_steps += reads_finished_in_this_chunk*global_context -> config.all_threads;

	unsigned long long guessed_all_steps = guessed_all_reads * (global_context -> index_block_number * 6 + 4);

	float finished_rate = finished_steps*1./guessed_all_steps;
	float reads_per_second = 0;

	if(task == STEP_VOTING)
		reads_per_second = finished_steps / (miltime() - global_context -> align_start_time) / (global_context -> index_block_number*6 + 4);
	else
		reads_per_second = finished_steps / (miltime() - global_context -> start_time) / (global_context -> index_block_number*6 + 4);
	//float exp_mins = (miltime() - global_context -> start_time) / finished_rate / 60; 

	//fprintf(stderr, "FINISHED=%llu, FINISHED_READS=%llu, ALL=%llu, ALLREADS=%llu, ALLCHUNKREADS=%llu; BEFORE_CHUK=%llu; CUR-BLK=%d; IND-BLK=%d\n", finished_steps, reads_finished_in_this_chunk, guessed_all_steps, guessed_all_reads,guessed_this_chunk_all_reads, guessed_all_reads_before_this_chunk,  global_context -> current_index_block_number , global_context -> index_block_number );

	if(current_read_no>1000 && !progress_report_callback)
		print_in_box(81,0,0, "%4d%%%% completed, %3d mins elapsed, total=%dk %s, rate=%2.1fk/s\r", (int)(finished_rate*100), (int)((miltime() - global_context -> start_time)/60),(int)(guessed_all_reads*1./1000), global_context -> input_reads.is_paired_end_reads?"frags":"reads", reads_per_second/1000, reads_finished_in_this_chunk);

	if(progress_report_callback)
	{
		progress_report_callback(10 ,task,(int)(10000*finished_rate));
		progress_report_callback(20 ,task,(int)(guessed_all_reads/1000));
	}
}


/*
int Xmain(int argc , char ** argv);

int main(int argc, char ** argv)
{
	int xk1=0;
	for(xk1=0; xk1<4; xk1++)
		Xmain(argc ,  argv);
	return 0;
}
*/


int parse_opts_core(int argc , char ** argv, global_context_t * global_context)
{
	int c;
	int option_index = 0;	

	optind = 1;
	opterr = 1;
	optopt = 63;

	while ((c = getopt_long (argc, argv, "ExsS:L:AHd:D:n:m:p:P:R:r:i:l:o:T:Q:I:t:B:b:Q:FcuUfM?", long_options, &option_index)) != -1)
	{
		switch(c)
		{
			case 'Q':
				global_context->config.multi_best_reads = atoi(optarg); 
				break;
			case 'H':
				global_context->config.use_hamming_distance_break_ties = 1;
				break;
			case 's':
				global_context->config.downscale_mapping_quality = 1;
				break;
			case 'M':
				global_context->config.do_big_margin_filtering_for_reads = 1;
				global_context->config.report_multi_mapping_reads = 0;
				break;

			case 'A':
				global_context->config.report_sam_file = 0;
				break;
			case 'E':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;

				break;
			case 'f':
				global_context->config.max_mismatch_exonic_reads = 200;
				global_context->config.max_mismatch_junction_reads = 200;
				global_context->config.do_fusion_detection = 1;
				global_context->config.minimum_subread_for_first_read = 1;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.total_subreads = 28;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 0;
				global_context->config.use_hamming_distance_in_exon = 1;
				break;
			case 'x':
				global_context->config.max_mismatch_exonic_reads = 10;
				global_context->config.max_mismatch_junction_reads = 1;
				global_context->config.ambiguous_mapping_tolerance = 39;
				global_context->config.extending_search_indels = 0;

				global_context->config.is_rna_seq_reads = 1;
				global_context->config.total_subreads = 14;
				global_context->config.minimum_subread_for_first_read = 3;
				global_context->config.minimum_subread_for_second_read = 1;
				global_context->config.high_quality_base_threshold = 990000;
				global_context->config.do_big_margin_filtering_for_junctions = 1;
				global_context->config.report_no_unpaired_reads = 0;
				global_context->config.limited_tree_scan = 1;
				global_context->config.use_hamming_distance_in_exon = 0;
				break;
			case 'S':
				global_context->config.is_first_read_reversed = optarg[0]=='r'?1:0;
				global_context->config.is_second_read_reversed = optarg[0]=='f'?0:1;
				break;
			case 'U':
				global_context->config.report_no_unpaired_reads = 1;
				break;
			case 'u':
				global_context->config.report_multi_mapping_reads = 0;
				break;
			case 'b':
				global_context->config.is_methylation_reads = 1;
				break;
			case 'D':
				global_context->config.maximum_pair_distance = atoi(optarg);
				break;
			case 'd':
				global_context->config.minimum_pair_distance = atoi(optarg);
				break;
			case 'n':
				global_context->config.total_subreads = atoi(optarg);
				break;
			case 'm':
				global_context->config.minimum_subread_for_first_read = atoi(optarg);
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

				if(global_context->config.max_indel_length <0)global_context->config.max_indel_length =0;
				if(global_context->config.max_indel_length > MAX_INSERTION_LENGTH)global_context->config.max_indel_length = MAX_INSERTION_LENGTH;
				if(global_context->config.max_indel_length > 16)
				{
					global_context->config.reassembly_subread_length = 12;
					global_context->config.reassembly_window_multiplex = 3;
					global_context->config.reassembly_start_read_number = 5;
					global_context->config.reassembly_tolerable_voting = 0;
					global_context->config.reassembly_window_alleles = 2;
					global_context->config.reassembly_key_length = 28;

					global_context->config.is_third_iteration_running = 1;
					global_context->config.max_mismatch_exonic_reads = 2;
					global_context->config.max_mismatch_junction_reads = 2;
					global_context->config.total_subreads = 28;
					global_context->config.do_big_margin_filtering_for_reads = 1;

					global_context->config.do_superlong_indel_detection = 0;
				}
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
			case 'B':
				global_context->config.is_first_iteration_running = 0;
				strcpy(global_context->config.medium_result_prefix, optarg);
				break;
			case 'c':
				global_context->config.space_type = GENE_SPACE_COLOR; 
				break;
				
			case 0:
				//if(strcmp("memory-optimisation",long_options[option_index].name)==0)
				//	memory_optimisation = atoi(optarg);
				break;
			case '?':
			default:
				return -1 ;
		}
	}


	return 0;
}

#ifdef MAKE_STANDALONE
int main_back(int argc , char ** argv)
{
	progress_report_callback = NULL;
#else
int subread_core_main(int argc , char ** argv)
{
	if(progress_report_callback) progress_report_callback(0,0,1);
#endif

	return core_main(argc, argv, parse_opts_core);
}


int check_configuration(global_context_t * global_context)
{
	int expected_type = FILE_TYPE_FAST_;
	if(global_context -> config.is_SAM_file_input && global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_BAM;
	else if(global_context -> config.is_SAM_file_input && !global_context -> config.is_BAM_input)
		expected_type = FILE_TYPE_SAM;
	else if(global_context -> config.is_gzip_fastq)
		expected_type = FILE_TYPE_GZIP_FAST_;
	
	if(global_context -> config.max_indel_length > 16)
		warning_file_limit();

	int wret = warning_file_type(global_context -> config.first_read_file, expected_type), wret2 = 0;
	if(global_context -> config.second_read_file[0])
	{
		if(expected_type==FILE_TYPE_FAST_ || expected_type==FILE_TYPE_GZIP_FAST_)
			wret2 = warning_file_type(global_context -> config.second_read_file, expected_type);
		else
			print_in_box(80,0,0,"Only one input SAM or BAM file is needed. The second input file is ignored.");
	}

	if(wret == -1 || wret2 == -1){
		return -1;
	}

	return 0;

}

int core_main(int argc , char ** argv, int (parse_opts (int , char **, global_context_t * )))
{
	//int memory_optimisation = 0;
	global_context_t * global_context;
	global_context = (global_context_t*)malloc(sizeof(global_context_t));
	init_global_context(global_context);


	int ret = parse_opts(argc , argv, global_context);
	if(ret) return ret;
	//global_context->config.reads_per_chunk = 200*1024;

	if(global_context->config.max_indel_length > 20 && !global_context->input_reads.is_paired_end_reads)
	{
		global_context->config.total_subreads = 28;
		global_context->config.reassembly_start_read_number = 3;
		global_context->config.do_superlong_indel_detection = 1;
	}


	ret = print_configuration(global_context);

	ret = ret || check_configuration(global_context);
	ret = ret || load_global_context(global_context);
	ret = ret || init_modules(global_context);


	ret = ret || read_chunk_circles(global_context);
	ret = ret || write_final_results(global_context);
	ret = ret || destroy_modules(global_context);
	ret = ret || destroy_global_context(global_context);
	ret = ret || show_summary(global_context);

	free(global_context);

	return ret;
}

// the new file name is written into fname then.

int convert_BAM_to_SAM(global_context_t * global_context, char * fname, int is_bam)
{
	char temp_file_name[200], *fline=malloc(3000), tmp_readname[MAX_READ_NAME_LEN+1];
	short tmp_flags;
	SamBam_FILE * sambam_reader;

	char * env_no_sort = getenv("SUBREAD_DO_NOT_CHECK_INPUT");
	if(env_no_sort){
		global_context->input_reads.is_paired_end_reads = 1;
		SUBREADprintf("\nWARNING: The SAM input file was assumed to be sorted by name.\nENV: SUBREAD_DO_NOT_CHECK_INPUT\n");
		return 0;
	}

	int is_file_sorted = 1;
	unsigned long long int read_no = 0;
	sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM:SAMBAM_FILE_SAM);
	if(!sambam_reader)
	{
		SUBREADprintf("ERROR: Unable to open file '%s'. The file is not in the specified format.\n", fname);
		free(fline);
		return -1;
	}
	while(1)
	{
		char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
		if(!is_ret) break;
		if(fline[0]=='@')continue;
		if(is_SAM_unsorted(fline, tmp_readname, &tmp_flags , read_no)){
			if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
			is_file_sorted = 0;
			read_no++;
			break;
		}
		read_no++;
		if(tmp_flags & 1) global_context->input_reads.is_paired_end_reads = 1;
		else global_context->input_reads.is_paired_end_reads = 0;

		if(! global_context->input_reads.is_paired_end_reads) break;
	}

	if(read_no<1)
	{
		SUBREADprintf("ERROR: unable to find any reads from file '%s'. The file is not in the specified format or is empty.\n", fname);
		return -1;
	}

	SamBam_fclose(sambam_reader);
	print_in_box(80,0,0,"The input %s file contains %s-end reads.", is_bam?"BAM":"SAM",  global_context->input_reads.is_paired_end_reads?"paired":"single");
	if(!is_file_sorted)
		print_in_box(80,0,0,"The input %s file is unsorted. Reorder it...", is_bam?"BAM":"SAM");
	else if(is_bam)
		print_in_box(80,0,0,"Convert the input BAM file...");

	if(is_bam || (global_context->input_reads.is_paired_end_reads && !is_file_sorted))
	{
		sprintf(temp_file_name, "%s.sam", global_context->config.temp_file_prefix);
		sambam_reader = SamBam_fopen(fname, is_bam?SAMBAM_FILE_BAM: SAMBAM_FILE_SAM);
		if(!sambam_reader){
			SUBREADprintf("Unable to open %s.\n", fname);
			return -1;
		}
		FILE * sam_fp = NULL;
		SAM_sort_writer writer;
		int writer_opened = 0;

		if(is_file_sorted) sam_fp = f_subr_open(temp_file_name,"w");
		else writer_opened = sort_SAM_create(&writer, temp_file_name, ".");

		if((is_file_sorted && !sam_fp) || (writer_opened && !is_file_sorted)){
			SUBREADprintf("Failed to write to the directory. You may not have permission to write to this directory or the disk is full.\n");
			return -1;
		}

		while(1)
		{
			char * is_ret = SamBam_fgets(sambam_reader, fline, 2999, 1);
			if(!is_ret) break;
			if(is_file_sorted)
				fputs(fline, sam_fp);
			else{
				int ret = sort_SAM_add_line(&writer, fline, strlen(fline));
				if(ret<0) {
					print_in_box(80,0,0,"ERROR: read name is too long; check input format.");
					break;
				}
			}
		}

		if(is_file_sorted) 
			fclose(sam_fp);
		else{
			sort_SAM_finalise(&writer);
			if(writer.unpaired_reads)
				print_in_box(80,0,0,"%llu single-end mapped reads in reordering.", writer.unpaired_reads);
		}

		SamBam_fclose(sambam_reader);
		strcpy(fname, temp_file_name);
		global_context -> will_remove_input_file = 1;
	}



	free(fline);
	return 0;
}

int convert_GZ_to_FQ(global_context_t * global_context, char * fname, int half_n)
{
	int is_OK = 0;
	char temp_file_name[200];
	char * linebuff=malloc(3001);
	gzFile * rawfp = gzopen(fname, "r");
	
	if(rawfp)
	{
		print_in_box(80,0,0,"Decompress %s...", fname);
		sprintf(temp_file_name, "%s-%d.fq", global_context->config.temp_file_prefix, half_n);
		FILE * outfp = fopen(temp_file_name, "w");
		if(outfp)
		{
			while(1)
			{
				char * bufr =gzgets(rawfp, linebuff, 3000);
				if(!bufr) break;
				fputs(bufr, outfp);
			}
			is_OK = 1;
			fclose(outfp);
		}
		else{
			SUBREADprintf("Unable to create temporary file '%s'\nThe program has to terminate.\nPlease run the program in a directory where you have the privilege to create files.\n", temp_file_name);
		}

		gzclose(rawfp);
	}

	strcpy(fname, temp_file_name);
	global_context -> will_remove_input_file |= (1<< (half_n-1));

	return !is_OK;
}

int core_geinput_open(global_context_t * global_context, gene_input_t * fp, int half_number, int is_init)
{
	char *fname;
	if(global_context->config.is_SAM_file_input)
	{

		fname = is_init?global_context ->config.first_read_file:global_context -> input_reads.first_read_file.filename;
		if(is_init && half_number == 1)
			if(convert_BAM_to_SAM(global_context, global_context ->config.first_read_file, global_context ->config.is_BAM_input)) return -1;
		if(!global_context->input_reads.is_paired_end_reads) half_number=0;
		return geinput_open_sam(fname, fp, half_number);

	}
	else
	{
		if(is_init)
		{
			if(global_context -> config.is_gzip_fastq)
				if(convert_GZ_to_FQ(global_context, (half_number==2)? global_context ->config.second_read_file : global_context ->config.first_read_file, half_number)) return -1;
			fname = (half_number == 2)?global_context -> config.second_read_file:global_context -> config.first_read_file;
		}
		else
			fname = (half_number == 2)?global_context -> input_reads.second_read_file.filename:global_context -> input_reads.first_read_file.filename;
		return geinput_open(fname, fp);
	}
}

void relocate_geinputs(global_context_t * global_context, thread_context_t * thread_context)
{
	if(thread_context)
	{
		thread_context -> reads_to_be_done = global_context -> input_reads.reads_in_blocks[thread_context -> thread_id];
		thread_context -> read_block_start = global_context -> input_reads.start_read_number_blocks[thread_context -> thread_id];

		thread_context -> ginp1 = (gene_input_t *)malloc(sizeof(gene_input_t));
		core_geinput_open(global_context, thread_context -> ginp1,1, 0);
		fseeko(thread_context -> ginp1 -> input_fp, global_context -> input_reads.first_file_blocks[thread_context -> thread_id], SEEK_SET);

		if(global_context -> input_reads.is_paired_end_reads)
		{
			thread_context -> ginp2 = (gene_input_t *)malloc(sizeof(gene_input_t));
			core_geinput_open(global_context, thread_context -> ginp2, 2, 0);
			fseeko(thread_context -> ginp2-> input_fp, global_context -> input_reads.second_file_blocks[thread_context -> thread_id], SEEK_SET);
		}
	}
}


int fetch_next_read_pair(global_context_t * global_context, thread_context_t * thread_context, gene_input_t* ginp1, gene_input_t* ginp2, int *read_len_1, int *read_len_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, char * qual_text_1, char *qual_text_2, int remove_color_head)
{
	int rl1, rl2=0;
	int is_second_R1, is_second_R2;

	do
	{
		is_second_R1 = 0; is_second_R2 = 0;
		rl1 = geinput_next_read_trim(ginp1, read_name_1, read_text_1 , qual_text_1, global_context->config.read_trim_5, global_context->config.read_trim_3, &is_second_R1);
		if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
		{
			if(isalpha(read_text_1[0]))
			{
				int xk1;
				for(xk1=2; read_text_1[xk1]; xk1++)
					read_text_1[xk1-2]=read_text_1[xk1];
				read_text_1[xk1-2]=0;
			}
		}

		if(ginp2)
		{
			rl2 = geinput_next_read_trim(ginp2, read_name_2, read_text_2 , qual_text_2, global_context->config.read_trim_5, global_context->config.read_trim_3, &is_second_R2);
			if(global_context->config.space_type == GENE_SPACE_COLOR && remove_color_head)
			{
				if(isalpha(read_text_2[0]))
				{
					int xk1;
					for(xk1=2; read_text_2[xk1]; xk1++)
						read_text_2[xk1-2]=read_text_2[xk1];
					read_text_2[xk1-2]=0;
				}
			}
		}
	} while(is_second_R1||is_second_R2) ;

	if( global_context->config.space_type == GENE_SPACE_COLOR)
	{
		rl1-=1;rl2-=1;
	}


	if(rl1>0 && (rl2>0 || !ginp2))
	{
		if(global_context->config.is_first_read_reversed)
		{
			reverse_read(read_text_1, rl1, global_context->config.space_type);
			if(qual_text_1)
				reverse_quality(qual_text_1, rl1);
		}

		if(ginp2 && global_context->config.is_second_read_reversed)
		{
			reverse_read(read_text_2, rl2, global_context->config.space_type);
			if(qual_text_2)
				reverse_quality(qual_text_2, rl2);
		}

		*read_len_1 = rl1;
		if(ginp2)
			*read_len_2 = rl2;
		return 0;
	}
	else	return 1;
}

int write_final_results(global_context_t * context)
{
	if(context -> config.output_prefix[0])
	{
		write_indel_final_results(context);

		if(context -> config.entry_program_name == CORE_PROGRAM_SUBJUNC)
			write_junction_final_results(context);

		if(context -> config.do_fusion_detection)
			write_fusion_final_results(context);
	}
	
	return 0;
}

int get_soft_clipping_length(char* CIGAR)
{
	int nch;
	int cigar_cursor, tmp_int = 0;

	for(cigar_cursor = 0; (nch = CIGAR[cigar_cursor]) > 0;cigar_cursor++)
	{
		if(isdigit(nch)) tmp_int = tmp_int*10+(nch-'0');
		else
		{
			if(nch=='S') return tmp_int;
			return 0;
		}
	}
	return 0;
}


#define CIGAR_PERFECT_SECTIONS 12

typedef struct{
	char current_cigar_decompress[CORE_MAX_CIGAR_STR_LEN + 1];
	char cigar [CORE_MAX_CIGAR_STR_LEN];

	unsigned int out_poses[CIGAR_PERFECT_SECTIONS];
	short out_lens[CIGAR_PERFECT_SECTIONS];
	char out_cigars[CIGAR_PERFECT_SECTIONS][60];
	char out_strands[CIGAR_PERFECT_SECTIONS];

	char additional_information[CORE_ADDITIONAL_INFO_LENGTH + 1];
	alignment_result_t * raw_result;

	unsigned int linear_position;
	unsigned int raw_linear;
	short soft_clipping_movements;
	char * chro;
	unsigned int offset;
	int strand;

	int mapping_quality;
}subread_output_tmp_t;

typedef struct{
	long long int *PE_distance;
	char * out_cigar_buffer[CIGAR_PERFECT_SECTIONS];
	subread_output_tmp_t *r1;
	subread_output_tmp_t *r2;
	subread_output_tmp_t ** out_pairs;
	alignment_result_t ** out_raws;

} subread_output_context_t;

void init_output_context(global_context_t * global_context ,subread_output_context_t * out_context)
{
	int xk1;
	out_context -> r1 = malloc(sizeof(subread_output_tmp_t) *global_context->config.multi_best_reads );

	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		out_context -> out_cigar_buffer[xk1] = malloc(60);

	out_context -> out_pairs = malloc(sizeof( subread_output_context_t *) * global_context->config.multi_best_reads * 2);
	out_context -> out_raws = malloc(sizeof( subread_output_context_t *) * global_context->config.multi_best_reads * 2);

	if(global_context -> input_reads.is_paired_end_reads)
	{
		out_context -> PE_distance = malloc(sizeof(long long) * global_context->config.multi_best_reads);
		out_context -> r2 = malloc(sizeof(subread_output_tmp_t) *global_context->config.multi_best_reads );
	}
	else{
		out_context -> PE_distance = NULL;
		out_context -> r2 = NULL;
	}
}

void destroy_output_context(global_context_t * global_context ,subread_output_context_t * out_context)
{
	int xk1;
	for(xk1=0;xk1<CIGAR_PERFECT_SECTIONS;xk1++)
		free(out_context -> out_cigar_buffer[xk1]);

	free(out_context -> out_raws);
	free(out_context -> out_pairs);
	free(out_context -> r1 );
	if(global_context -> input_reads.is_paired_end_reads)
	{
		free(out_context -> r2);
		free(out_context -> PE_distance);
	}
}

int locate_current_value_index(global_context_t * global_context, thread_context_t * thread_context, alignment_result_t * result, int rlen);
int calc_edit_dist(global_context_t * global_context, alignment_result_t * current_result, char * cigar , unsigned int pos ,  char * read_text)
{
	locate_current_value_index(global_context, NULL, current_result, 1);
	gene_value_index_t * current_value_index = global_context->current_value_index; 

	int cigar_cursor=0;
	unsigned int chro_cursor = pos, tmpi=0;
	int read_cursor = 0;
	int all_mm = 0;

	while(1)
	{
		char nch = cigar[cigar_cursor++];
		if(!nch)break;

		if('0'<=nch && '9' >=nch)
		{
			tmpi = tmpi*10+nch-'0';
		}
		else{
			if(nch == 'M')
			{
				int matched = match_chro(read_text + read_cursor, current_value_index, chro_cursor, tmpi, 0, global_context -> config.space_type);
				all_mm += tmpi - matched;
				chro_cursor += tmpi;
				read_cursor += tmpi;
			}
			else if(nch == 'N' || nch == 'D')
			{
				if('D' == nch) all_mm+=tmpi;
				chro_cursor += tmpi;
			}
			else if(nch == 'I')
			{
				read_cursor += tmpi;
				all_mm+=tmpi;
			}
			else if(nch == 'S')
			{
				chro_cursor += tmpi;
				read_cursor += tmpi;
			}
	
			tmpi = 0;
		}
	}
	return all_mm;
}

int convert_read_to_tmp(global_context_t * global_context , subread_output_context_t * output_context, int read_number, int is_second_read, int read_len, char * read_text, char * qual_text, alignment_result_t * current_result, subread_output_tmp_t * r)
{

	int is_r_OK;
	r -> raw_result = current_result;
	r -> additional_information[0]=0;

	is_r_OK = (current_result -> result_flags & CORE_IS_FULLY_EXPLAINED) > 0;

	if(is_r_OK)
		if((current_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
			is_r_OK = 0;


	if(is_r_OK){

		int current_strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
		int is_first_section_jumped = 0;


		if( current_result -> cigar_string[0] == -1)
		{
			bincigar2cigar( r-> current_cigar_decompress, CORE_MAX_CIGAR_STR_LEN, current_result -> cigar_string + 1, CORE_MAX_CIGAR_LEN, read_len);
			is_first_section_jumped = 1;
		}
		else
			bincigar2cigar( r-> current_cigar_decompress, CORE_MAX_CIGAR_STR_LEN, current_result -> cigar_string, CORE_MAX_CIGAR_LEN, read_len);

		int chimeric_sections = 0;
		int current_repeated_times;

		current_repeated_times = is_ambiguous_voting(global_context, read_number, is_second_read, current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, read_len, current_strand);

		r->raw_linear = current_result -> selected_position;
		r->linear_position = current_result -> selected_position;
		r->mapping_quality = current_result -> final_quality;
		if(current_repeated_times>1) r->mapping_quality = 0;

		strcpy(r->cigar, r -> current_cigar_decompress);
		r->strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

		//sprintf(r->additional_information, "\tSM:i:%d",current_result -> final_mismatched_bases);


		//printf("SM='%s'\n", r->additional_information);

		//#warning "THIS LINE CAN WORK VERY BADLY!!! DO NOT INCLUDE IT IN THE RELEASE UNTIL IT IS FULLY TESTED!"
		if(global_context -> config.SAM_extra_columns)
			sprintf(r->additional_information + strlen( r->additional_information), "\tSB:i:%d\tSC:i:%d\tSD:i:%d\tSN:i:%u\tSP:Z:%s", current_result -> used_subreads_in_vote, current_result -> selected_votes, current_result -> noninformative_subreads_in_vote, read_number, (current_result -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 

		if(global_context -> config.do_fusion_detection)
		{			
			chimeric_sections = chimeric_cigar_parts(global_context, r->linear_position, is_first_section_jumped ^ current_strand, is_first_section_jumped ,  r->current_cigar_decompress , r->out_poses, output_context->out_cigar_buffer, r->out_strands, read_len, r->out_lens);

			int xk1;
			strcpy(r->out_cigars[0], output_context->out_cigar_buffer[0]);
			for(xk1=1; xk1<chimeric_sections; xk1++)
			{
				unsigned int chimeric_pos;
				char * chimaric_chr;
				strcpy(r->out_cigars[xk1], output_context->out_cigar_buffer[xk1]);

				if(0==locate_gene_position_max(r->out_poses[xk1],& global_context -> chromosome_table, & chimaric_chr, & chimeric_pos, 0+r->out_lens[xk1]))
				{
					int soft_clipping_movement = 0;
					soft_clipping_movement = get_soft_clipping_length(r->out_cigars[xk1]);
					char strand_xor = (r->out_strands[xk1] == '-');
					assert(chimaric_chr);
					sprintf(r->additional_information + strlen(r->additional_information), "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s", r->out_cigars[xk1] , chimeric_pos + soft_clipping_movement + 1, strand_xor?'-':'+' , chimaric_chr );
				}
				else is_r_OK = 0;
			}
			r->linear_position = r->out_poses[0];
			r->strand = r->out_strands[0]=='-';

			strcpy(r->cigar , r->out_cigars[0]);
		}
		r->soft_clipping_movements = get_soft_clipping_length(r->cigar);
	}
	else if(global_context -> config.SAM_extra_columns)
	{
		sprintf(r->additional_information + strlen( r->additional_information), "\tSB:i:%d\tSC:i:%d\tSD:i:%d\tSN:i:%u\tSP:Z:%s", current_result -> used_subreads_in_vote, current_result -> selected_votes, current_result -> noninformative_subreads_in_vote, read_number, (current_result -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 
	}

	if(is_r_OK)
	{
		if(locate_gene_position_max(r->linear_position,& global_context -> chromosome_table, &r-> chro , &r -> offset, read_len))
		{
			is_r_OK = 0;
		}
		else
		{
			r -> offset++;
			assert(r-> chro);
		}

		if(global_context -> config.is_rna_seq_reads && !(current_result -> result_flags & CORE_NOTFOUND_DONORS))
		{
			sprintf(r->additional_information + strlen(r->additional_information), "\tXS:A:%c", (current_result -> result_flags & CORE_IS_GT_AG_DONORS)?'+':'-');
		}

		/*
		if(global_context -> config.more_accurate_fusions)
		{
			
			sprintf(r->additional_information + strlen(r->additional_information), "\tDF:i:%d", current_result -> best_second_diff_bases >0?current_result -> best_second_diff_bases:9999 );
		}*/


	}

	return is_r_OK;
}


int calculate_fragment_combinations(global_context_t * global_context, subread_output_context_t * out_context, int read_number , int read_len_1, int read_len_2, char * read_name_1, char * read_name_2, char * read_text_1, char * read_text_2, char * qual_text_1, char * qual_text_2)
{
	int mapping_location_number;

	// Features that need attention:

	// (1) fusion cigar split
	// (2) soft-clipping moving mapping location
	// (3) making the two ends to have the best strand: if both ends mapped, as what they should be; 
	//		       if only one end is mapped: the other end is flipped to the same strand.
	// (4) ambiguous voting

	// Additional criteria: 
	// (1) a read must be entirely in a chromosome, or it is reported as unmapped.

	// The return value is the number of locations.
	// It must be greater than 0, even both reads were unmapped.
	int read_1_locations = 0;
	int read_2_locations = 0;

	if(global_context -> config.report_multiple_best_in_pairs)
	{
		int retX = 0, is_two_read_mapped = 0;
		int is_2_OK = 0, is_1_OK = 0;


		for( mapping_location_number = 0 ; mapping_location_number < global_context->config.multi_best_reads; mapping_location_number++)
		{
			assert(global_context->input_reads.is_paired_end_reads);	// this function is not provided for the public.
			alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0,  mapping_location_number);
			is_1_OK = convert_read_to_tmp(global_context, out_context, read_number , 0 , read_len_1, read_text_1, qual_text_1, current_result, &out_context->r1[read_1_locations]);

			alignment_result_t * current_result2 = _global_retrieve_alignment_ptr(global_context  , read_number, 1,  mapping_location_number);
			is_2_OK = convert_read_to_tmp(global_context , out_context, read_number , 1, read_len_2, read_text_2, qual_text_2, current_result2, &out_context->r2[read_2_locations]);

			
			if(mapping_location_number == 0)
				is_two_read_mapped = is_1_OK && is_2_OK;

			if(mapping_location_number == 0 || (is_1_OK && is_2_OK))
				retX++;
			else break;

			if(is_1_OK)
			{
				out_context -> out_pairs[mapping_location_number*2  ]  = out_context->r1 + read_1_locations;
				assert(out_context -> out_pairs[mapping_location_number*2  ]->chro);
			}
			if(is_2_OK)
			{
				out_context -> out_pairs[mapping_location_number*2+1]  = out_context->r2 + read_2_locations;
				assert(out_context -> out_pairs[mapping_location_number*2+1]->chro);
			}

			if(global_context -> config.SAM_extra_columns)
			{
				out_context -> out_raws[mapping_location_number*2] = current_result;
				out_context -> out_raws[mapping_location_number*2+1] = current_result2;
			}

			if(is_1_OK)read_1_locations++;
			if(is_2_OK)read_2_locations++;

			if(!is_two_read_mapped) break;
		}


		return retX;
	}

	// First, counting the number of mapping locations for each end.
	// 
	for( mapping_location_number = 0 ; mapping_location_number < global_context->config.multi_best_reads; mapping_location_number++)
	{
		alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0,  mapping_location_number);

		int is_2_OK = 0, is_1_OK = 0, xx1, is_fresh;

		is_fresh = 1;
		for(xx1=0;xx1 < read_1_locations; xx1++)
		{

			if(current_result -> selected_position == out_context->r1[xx1].raw_linear)
				is_fresh = 0;
		}

		//printf("SB:i:%d\tSC:i:%d\tSD:i:%d\tSP:Z:%s\n", current_result -> used_subreads_in_vote, current_result -> selected_votes, current_result -> noninformative_subreads_in_vote, (current_result -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 
		if(is_fresh)
			is_1_OK = convert_read_to_tmp(global_context, out_context, read_number , 0 , read_len_1, read_text_1, qual_text_1, current_result, &out_context->r1[read_1_locations]);

		if(global_context->input_reads.is_paired_end_reads)
		{
			alignment_result_t * current_result2 = _global_retrieve_alignment_ptr(global_context  , read_number, 1,  mapping_location_number);


			is_fresh = 1;
			for(xx1=0;xx1 < read_2_locations; xx1++)
				if(current_result2 -> selected_position == out_context->r2[xx1].linear_position)
					is_fresh = 0;

			if(is_fresh)
				is_2_OK = convert_read_to_tmp(global_context , out_context, read_number , 1, read_len_2, read_text_2, qual_text_2, current_result2, &out_context->r2[read_2_locations]);

			if(is_2_OK)read_2_locations++;
		}
		//printf("L1,2_OK = %d, %d\n", is_1_OK, is_2_OK);

		if(is_1_OK)read_1_locations++;

		if(! ( is_2_OK || is_1_OK))
			break;
	}


	// now all the potential locations were written into out_context->r1 and out_context->r2.
	// we find the best combinations between out_context->r1 and out_context->r2 now.
	int ret = 0;
	if(global_context-> config.multi_best_reads == 1)
	{
		out_context -> out_pairs[0]   =  read_1_locations?out_context -> r1 : NULL;
		out_context -> out_pairs[1]   =  read_2_locations?out_context -> r2 : NULL;
		ret = 1;
	}
	else if(global_context->input_reads.is_paired_end_reads)
	{
		int r1_xx, r2_xx, yy, zz;
		memset(out_context->PE_distance, 0x7F, sizeof(long long) * global_context->config.multi_best_reads);


		out_context -> out_pairs[0] = NULL;
		out_context -> out_pairs[1] = NULL;

		if(read_1_locations > 0 && read_2_locations > 0)
		{
			for(r1_xx=0; r1_xx < read_1_locations; r1_xx ++)
				for(r2_xx=0; r2_xx < read_2_locations; r2_xx ++)
				{
					long long int distance = 0x100000000ll;

					// the chro in the data structure is a pointer to the offset table; they can be compaired directly.
					if(out_context->r1 [r1_xx] . chro == out_context->r2 [r2_xx].chro)
					{
						distance = out_context->r1 [r1_xx] . linear_position;
						distance -= out_context->r2 [r2_xx] . linear_position;
						if(0>distance) distance=-distance;
					}

					if(distance < out_context -> PE_distance[global_context->config.multi_best_reads - 1])
					{
						for(yy = 0; yy < global_context->config.multi_best_reads; yy++)
						{
							if(distance < out_context -> PE_distance[yy])
								break;
						}

						// now replace the 2*yy-th and (2*yy+1)-th items in array output_records with the current couple.

						for(zz = global_context->config.multi_best_reads - 2 ; zz >= yy; zz--)
						{
							out_context -> PE_distance[zz+1] = out_context -> PE_distance[zz];
							memcpy(&(out_context -> out_pairs[zz*2 + 2]), &(out_context -> out_pairs[zz*2]), 2*sizeof(void *));
						}
						if(yy < global_context->config.multi_best_reads)
						{
							out_context -> PE_distance[yy] = distance;
							out_context -> out_pairs[yy*2]   = out_context -> r1  + r1_xx;
							out_context -> out_pairs[yy*2+1] = out_context -> r2  + r2_xx;

							ret ++;
						}
						if(ret > global_context->config.multi_best_reads) ret = global_context->config.multi_best_reads;
					}
				}
		}
		else
		{
			for(yy = 0; yy < max(read_1_locations, read_2_locations); yy++)
			{
				out_context -> out_pairs[yy*2  ] = read_1_locations == 0? NULL :  out_context -> r1  + yy;
				out_context -> out_pairs[yy*2+1] = read_2_locations == 0? NULL :  out_context -> r2  + yy;
			}
			ret =  max(read_1_locations, read_2_locations);
		}

//		if(memcmp("V0112_0155:7:1102:10778:2461", read_name_1, 27)==0)
//			SUBREADprintf("NNNLLL=%d,%d,%d\n", read_1_locations, read_2_locations, ret);
	}
	else
	{
		int yy;
		ret = read_1_locations;
		for(yy = 0; yy < read_1_locations; yy++)
			out_context -> out_pairs[yy  ] =  out_context -> r1 +yy;
	}
	//assert(ret <=global_context->config.multi_best_reads);
	return max(1, ret);
}

int calc_tlen(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int read_len_1, int  read_len_2);

int calc_flags(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int is_second_read, int read_len_1, int read_len_2, int current_location_no)
{
	int ret;

	if(global_context->input_reads.is_paired_end_reads)
	{
		ret  = SAM_FLAG_PAIRED_TASK;
		ret |= is_second_read?SAM_FLAG_SECOND_READ_IN_PAIR:SAM_FLAG_FIRST_READ_IN_PAIR;

		subread_output_tmp_t * this_rec = is_second_read?rec2:rec1;
		subread_output_tmp_t * mate_rec = is_second_read?rec1:rec2;

		if(this_rec == NULL)  ret |= SAM_FLAG_UNMAPPED;
		else
			if(this_rec->strand + is_second_read == 1) ret |= SAM_FLAG_REVERSE_STRAND_MATCHED;

		if(mate_rec == NULL)  ret |= SAM_FLAG_MATE_UNMATCHED;
		else
			if(mate_rec->strand + is_second_read != 1) ret |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;

		if(rec1 && rec2)
		{
			int TLEN = calc_tlen(global_context , rec1, rec2, read_len_1, read_len_2);
			int is_PEM = 0;
			if(TLEN >= global_context->config. minimum_pair_distance  && TLEN <= global_context-> config.maximum_pair_distance &&  this_rec->strand == mate_rec->strand)
			{
				if(global_context -> config.is_first_read_reversed && !(global_context -> config.is_second_read_reversed))
				{
					if(this_rec -> strand == 0)
					{
						if((is_second_read + (mate_rec-> offset > this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset)
							is_PEM = 1;
					}
				}
				else
				{
					if(this_rec -> strand)
					{
						if((is_second_read + (mate_rec-> offset < this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset) is_PEM = 1;
					}else
					{
						if((is_second_read + (mate_rec-> offset > this_rec -> offset) == 1) || mate_rec-> offset == this_rec -> offset) is_PEM = 1;

					}
				}
			}
			if(is_PEM) ret |= SAM_FLAG_MATCHED_IN_PAIR;
		}
	}
	else
	{
		ret = 0;

		if(rec1 == NULL)  ret |= SAM_FLAG_UNMAPPED;
		else
			if(rec1->strand) ret |= SAM_FLAG_REVERSE_STRAND_MATCHED;
	}

	if(current_location_no>0){
		if((rec1 && !is_second_read)||(rec2 && is_second_read))
			ret |= SAM_FLAG_SECONDARY_MAPPING;
	}

	return ret;
}


int calc_tlen(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int read_len_1, int  read_len_2)
{
	long long int ret = rec1->offset;
	ret -= rec2->offset;
	if(ret<0)ret=-ret;

	if(rec1->offset - rec1->soft_clipping_movements < rec2->offset - rec2->soft_clipping_movements ) ret += read_len_2;
	else ret += read_len_1;

	ret = max(read_len_1, ret);
	ret = max(read_len_2, ret);
	return (int)ret;
}

int calc_should_reverse(global_context_t * global_context, subread_output_tmp_t * rec1 , subread_output_tmp_t * rec2, int is_second_read)
{
	subread_output_tmp_t * cur_rec = is_second_read?rec2:rec1;

	if(cur_rec)
		return cur_rec->strand;
	else
	{
		return is_second_read;
		//subread_output_tmp_t * mate_rec = is_second_read?rec2:rec1;
		//return mate_rec?mate_rec->strand:is_second_read;
	}
}

void remove_nm_i(char * st)
{
	char * st_nm = strstr(st, "\tNM:i:");
	if(st_nm)
	{
		char * write_cursor = st_nm;
		int started = 0;
		st_nm++;
		while(1)
		{
			char nch = *st_nm;
			if(!nch)
			{
				*write_cursor=0;
				break;
			}
			if(nch == '\t') started = 1;
			if(started)
			{
				*write_cursor=*st_nm;
				write_cursor++;
			}
			st_nm++;
		}
	}
}

#define write_chunk_results_145 write_chunk_results

// rec1 or rec2 is OK if they are not NULL.
void write_single_fragment(global_context_t * global_context, subread_output_tmp_t * rec1, alignment_result_t * raw_r1, subread_output_tmp_t * rec2, alignment_result_t * raw_r2, int all_locations , int current_location , char * read_name_1, char * read_name_2, int read_len_1, int read_len_2, char * read_text_1, char * read_text_2, char * qual_text_1, char * qual_text_2, int * is_read1_reversed, int * is_read2_reversed )
{

//	if( all_locations < 2) return;


	int flag1 = calc_flags( global_context , rec1, rec2, 0,  read_len_1,  read_len_2, current_location);

	int flag2 = -1;

	if(global_context->input_reads.is_paired_end_reads)
	{
		flag2 = calc_flags( global_context , rec1, rec2, 1,  read_len_1,  read_len_2, current_location);
		if((0 == current_location) && (flag2 & SAM_FLAG_MATCHED_IN_PAIR)) global_context->all_correct_PE_reads  ++;
	}

	int tlen = 0;

	// rec -> chro is a pointer to the offset table; the pointers can be compared.
	if(rec1 && rec2 && rec1->chro == rec2->chro)tlen = calc_tlen(global_context , rec1, rec2,  read_len_1,  read_len_2);



	int applied_reverse_space;
	applied_reverse_space = global_context->config.space_type;
	if(global_context -> config.convert_color_to_base)
	{
		colorread2base(read_text_1, read_len_1+1);
		if(global_context->input_reads.is_paired_end_reads)
			colorread2base(read_text_2, read_len_2+1);
		applied_reverse_space = GENE_SPACE_BASE;
	}

	int should_1_reverse = calc_should_reverse( global_context , rec1, rec2, 0);

	if(should_1_reverse + (*is_read1_reversed) == 1)
	{
		reverse_read(read_text_1, read_len_1 + global_context->config.convert_color_to_base, applied_reverse_space);
		reverse_quality(qual_text_1, read_len_1);
		(*is_read1_reversed) = !(*is_read1_reversed);
	}

	if(global_context->input_reads.is_paired_end_reads)
	{
		int should_2_reverse = calc_should_reverse( global_context , rec1, rec2, 1);

		if(should_2_reverse + (*is_read2_reversed) == 1)
		{
			reverse_read(read_text_2, read_len_2 + global_context->config.convert_color_to_base, applied_reverse_space);
			reverse_quality(qual_text_2, read_len_2);
			(*is_read2_reversed) = !(*is_read2_reversed);
		}
	}
	remove_backslash(read_name_1);
	remove_backslash(read_name_2);

	int display_offset1 = 0, display_tailgate1 = 0;
	int display_offset2 = 0, display_tailgate2 = 0;

	if(global_context -> config.space_type == GENE_SPACE_COLOR)
	{

		if(rec1 && rec1 -> strand)
		{
			display_offset1 = 0;
			display_tailgate1 = 1;
		}
		else
		{
			if(!global_context -> config.convert_color_to_base)
			{
				// the first base was a fake prime base; the second base is the first meaningful base.
				char second_char = read_text_1[1];
				read_text_1[1] = color2char(second_char, read_text_1[0]);
			}
			display_offset1 = 1;
			display_tailgate1 = 0;
		}

		if(rec2 && rec2 -> strand)
		{
			if(!global_context -> config.convert_color_to_base)
			{
				// the first base was a fake prime base; the second base is the first meaningful base.
				char second_char = read_text_2[1];
				read_text_2[1] = color2char(second_char, read_text_2[0]);
			}
			display_offset2 = 1;
			display_tailgate2 = 0;
		}
		else
		{
			display_offset2 = 0;
			display_tailgate2 = 1;
		}
		assert(display_offset1 + display_tailgate1 == 1);
		assert(display_offset2 + display_tailgate2 == 1);
	}
	if(!qual_text_1[0])
	{
		int xi2;
		for(xi2=display_offset1;read_text_1[xi2 + display_tailgate1];xi2++) qual_text_1[xi2 - display_offset1] = 'I';
		qual_text_1[xi2 - display_offset1]=0;
		for(xi2=display_offset2;read_text_2[xi2 + display_tailgate2];xi2++) qual_text_2[xi2 - display_offset2] = 'I';
		qual_text_2[xi2 - display_offset2]=0;

	}



	if(display_tailgate1)
		read_text_1[strlen(read_text_1)-1]=0;
	if(display_tailgate2)
		read_text_2[strlen(read_text_2)-1]=0;


	if(global_context -> config.space_type == GENE_SPACE_BASE){ 
		if(rec1 && !strstr(rec1->additional_information, "\tNM:i:")){
			short rec1_edit = calc_edit_dist(global_context, rec1->raw_result, rec1->cigar , rec1->linear_position, read_text_1);
			sprintf(rec1->additional_information + strlen( rec1->additional_information), "\tNM:i:%d", rec1_edit );
		}
		if(global_context->input_reads.is_paired_end_reads && rec2 && !strstr(rec2->additional_information, "\tNM:i:"))
		{
			short rec2_edit = calc_edit_dist(global_context, rec2->raw_result, rec2->cigar , rec2->linear_position, read_text_2);
			sprintf(rec2->additional_information + strlen( rec2->additional_information), "\tNM:i:%d", rec2_edit );
		}
	}

	char extra_additional_1 [200+CORE_ADDITIONAL_INFO_LENGTH], extra_additional_2[200+CORE_ADDITIONAL_INFO_LENGTH];

	extra_additional_1[0]=0;
	extra_additional_2[0]=0;

	if(global_context -> config.SAM_extra_columns)
	{
		if(raw_r1)
			sprintf(extra_additional_1, "SB:i:%d\tSC:i:%d\tSD:i:%d\tSP:Z:%s\t", raw_r1 -> used_subreads_in_vote, raw_r1 -> selected_votes, raw_r1 -> noninformative_subreads_in_vote, (raw_r1 -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 
		if(raw_r2)
			sprintf(extra_additional_2, "SB:i:%d\tSC:i:%d\tSD:i:%d\tSP:Z:%s\t", raw_r2 -> used_subreads_in_vote, raw_r2 -> selected_votes, raw_r2 -> noninformative_subreads_in_vote, (raw_r2 -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 
	}


	if(rec1)
		sprintf(extra_additional_1 +strlen(extra_additional_1), "HI:i:%d\tNH:i:%d", current_location+1, all_locations);
	if(rec2)
		sprintf(extra_additional_2 +strlen(extra_additional_2), "HI:i:%d\tNH:i:%d", current_location+1, all_locations);

	if(global_context->config.read_group_id[0])
	{
		snprintf(extra_additional_1+strlen(extra_additional_1),100, "\tRG:Z:%s", global_context->config.read_group_id);
		snprintf(extra_additional_2+strlen(extra_additional_2),100, "\tRG:Z:%s", global_context->config.read_group_id);
	}

	char * out_chro1, * out_chro2, *out_cigar1, *out_cigar2;
	if(rec1)
	{
		assert(rec1->chro);
		strcat(extra_additional_1, rec1->additional_information);
		assert(rec1->chro);
		//printf("STRCAT: + '%s' = '%s'\n",  rec1->additional_information, extra_additional_1);
		out_chro1 = rec1->chro;
		out_cigar1 = rec1->cigar;
	}
	else
	{
		out_chro1 = "*";
		out_cigar1 = "*";
	}

	if(rec2)
	{
		assert(rec2->chro);
		strcat(extra_additional_2, rec2->additional_information);
		assert(rec2->chro);
		out_chro2 = rec2->chro;
		out_cigar2 = rec2->cigar;
	}
	else
	{
		out_chro2 = "*";
		out_cigar2 = "*";
	}

	int out_offset1=0, out_offset2=0;
	long long int out_tlen1, out_tlen2;

	out_tlen1 = tlen;
	out_tlen2 = tlen;
	if(rec1 && rec2)
	{
		if( rec1->offset >  rec2->offset) out_tlen1 = - out_tlen1;
		else	out_tlen2 = -out_tlen2;
	}

	if(0==current_location)
	{
		if(global_context -> input_reads.is_paired_end_reads)
		{
			if(rec1 || rec2)
				global_context -> all_mapped_reads++;
		}
		else if(rec1)
			global_context -> all_mapped_reads++;
	}

	if(rec1)
		out_offset1 = rec1->offset + rec1 -> soft_clipping_movements;
	if(rec2)
		out_offset2 = rec2->offset + rec2 -> soft_clipping_movements;


	int  out_mapping_quality1 = 0, out_mapping_quality2 = 0;
	if(rec1)
		 out_mapping_quality1 = rec1->mapping_quality;
	if(rec2)
		 out_mapping_quality2 = rec2->mapping_quality;

	char * mate_chro_for_1  = out_chro2;
	char * mate_chro_for_2  = out_chro1;

	if(out_chro1 == out_chro2 && out_chro1 && out_chro1[0]!='*'){
		mate_chro_for_1="=";
		mate_chro_for_2="=";
	}


	if(global_context -> config.is_BAM_output){
		SamBam_writer_add_read(global_context -> output_bam_writer, read_name_1, flag1,  out_chro1 , out_offset1, out_mapping_quality1, out_cigar1, out_chro2 , out_offset2, out_tlen1, read_len_1, read_text_1 + display_offset1, qual_text_1, extra_additional_1);

		if(global_context->input_reads.is_paired_end_reads)
			SamBam_writer_add_read(global_context -> output_bam_writer, read_name_2, flag2,  out_chro2 , out_offset2, out_mapping_quality2, out_cigar2, out_chro1 , out_offset1, out_tlen2, read_len_2, read_text_2 + display_offset2, qual_text_2, extra_additional_2);
	}
	else
	{
		sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s%s\n", read_name_1, flag1, out_chro1, out_offset1, out_mapping_quality1, out_cigar1, mate_chro_for_1, out_offset2, out_tlen1, read_text_1 + display_offset1, qual_text_1, extra_additional_1[0]?"\t":"", extra_additional_1);
		if(global_context->input_reads.is_paired_end_reads)
			sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s%s\n", read_name_2, flag2, out_chro2, out_offset2, out_mapping_quality2, out_cigar2, mate_chro_for_2, out_offset1, out_tlen2, read_text_2 + display_offset2, qual_text_2, extra_additional_2[0]?"\t":"", extra_additional_2);
	}
}

int write_chunk_results_145(global_context_t * global_context)
{
	
	/**********************************************
	 **********************************************
	 ** Initiate the memory blocks
	 **********************************************
	 **********************************************
	 */
	gene_input_t * ginp1, * ginp2=NULL;

	ginp1 = &global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		ginp2 = &global_context->input_reads.second_read_file;

	int read_number, sqr_read_number = 0, sqr_interval = global_context -> processed_reads_in_chunk/10; 
	
	subread_output_context_t out_context;
	init_output_context(global_context, &out_context);	

	char * read_text_1, * read_text_2;
	char * qual_text_1, * qual_text_2;
	char * read_name_1, * read_name_2;
	char * output_line_buffer;

	int read_len_1 = 0, read_len_2 = 0;
	read_text_1 = malloc(sizeof(char) * (MAX_READ_LENGTH+1));
	read_text_2 = malloc(sizeof(char) * (MAX_READ_LENGTH+1));
	qual_text_1 = malloc(sizeof(char) * (MAX_READ_LENGTH+1));
	qual_text_2 = malloc(sizeof(char) * (MAX_READ_LENGTH+1));
	read_name_1 = malloc(sizeof(char) * (MAX_READ_NAME_LEN+1));
	read_name_2 = malloc(sizeof(char) * (MAX_READ_NAME_LEN+1));
	output_line_buffer = malloc(sizeof(char) * ( 2* MAX_READ_LENGTH + 2 * MAX_CHROMOSOME_NAME_LEN + MAX_READ_NAME_LEN + CORE_MAX_CIGAR_STR_LEN + CORE_ADDITIONAL_INFO_LENGTH + 100));

	/**********************************************
	 **********************************************
	 ** Going through all reads in this chunk.
	 **********************************************
	 **********************************************
	 */
	for(read_number = 0; read_number < global_context -> processed_reads_in_chunk ; read_number++)
	{		
		int output_alignment_number, is_read1_reversed, is_read2_reversed;

		if(sqr_read_number > sqr_interval)
		{
			show_progress(global_context, NULL, read_number, STEP_WRITE_CHUNK_RESULTS);
			sqr_read_number = 0;
		}

		sqr_read_number ++;
		fetch_next_read_pair(global_context, NULL, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 0);
		if(qual_text_1[0] && global_context -> config.phred_score_format == FASTQ_PHRED64)
		{
			fastq_64_to_33(qual_text_1);
			if(global_context->input_reads.is_paired_end_reads)
				fastq_64_to_33(qual_text_2);
		
		}

		is_read1_reversed = 0;
		is_read2_reversed = 0;

		memset( out_context.out_pairs, 0, sizeof(void *) *2* global_context->config.multi_best_reads);
		memset( out_context.out_raws, 0, sizeof(void *) *2* global_context->config.multi_best_reads);

		// Array output_records is organised in the order of [R1, R2, R1, R2, R1, R2...]. The total number of items is output_fragment_combinations (*2 for paired-end reads).
		int output_fragment_combinations = calculate_fragment_combinations(global_context , &out_context, read_number, read_len_1, read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2);


		/**********************************************
		 **********************************************
		 ** Write the fragment results. Note that each fragment MUST be paired.
		 **********************************************
		 **********************************************
		 */

	//	if(output_fragment_combinations >= 0) printf("MM READ=%d : %s\n",  output_fragment_combinations, read_name_1); 

		for(output_alignment_number = 0; output_alignment_number < max(1,output_fragment_combinations);output_alignment_number ++)
		{
			char * read_text_1_raw =NULL, * read_text_2_raw =NULL, is_null_qual = 0;
			if(global_context -> config.space_type == GENE_SPACE_COLOR)
			{
				read_text_1_raw = malloc(1250);
				read_text_2_raw = malloc(1250);
				is_null_qual = qual_text_1[0]==0;
				strcpy(read_text_1_raw, read_text_1);
				strcpy(read_text_2_raw, read_text_2);
				// for color-space reads, the read text is re-written everytime because it is impossible to re-reverse the reads.
			}

	

			
			subread_output_tmp_t * rec1 = out_context.out_pairs[2*output_alignment_number];
			subread_output_tmp_t * rec2 = out_context.out_pairs[2*output_alignment_number+1];
			
			
			if(rec1)
			{
			//	SUBREADprintf("OUTN=%d; selected_votes=%d\n\n", output_alignment_number, rec1->raw_result -> selected_votes);
				assert(rec1->chro);
			}
			if(rec2)
			{
				assert(rec2->chro);
			//	SUBREADprintf("OUTT=%d; selected_votes=%d\n\n", output_alignment_number, rec2->raw_result -> selected_votes);
			}

			alignment_result_t * raw_r1, * raw_r2;

			if(rec1 || rec2  || !global_context->config.ignore_unmapped_reads)
			{
				if(global_context->input_reads.is_paired_end_reads)
				{
					raw_r1 = out_context.out_raws[2*output_alignment_number];
					raw_r2 = out_context.out_raws[2*output_alignment_number+1];
					write_single_fragment(global_context, out_context.out_pairs[2*output_alignment_number], raw_r1, out_context.out_pairs[2*output_alignment_number + 1], raw_r2, output_fragment_combinations , output_alignment_number , read_name_1, read_name_2, read_len_1, read_len_2, read_text_1, read_text_2, qual_text_1, qual_text_2, &is_read1_reversed, &is_read2_reversed );
				}
				else
				{
					raw_r1 = out_context.out_raws[output_alignment_number];
					write_single_fragment(global_context, out_context.out_pairs[output_alignment_number], raw_r1, NULL, NULL, output_fragment_combinations , output_alignment_number , read_name_1, read_name_2, read_len_1, read_len_2, read_text_1, read_text_2, qual_text_1, qual_text_2, &is_read1_reversed, &is_read2_reversed );
				}
			}

			if(global_context -> config.space_type == GENE_SPACE_COLOR)
			{
				// for color-space reads, the read text is re-written everytime because it is impossible to re-reverse the reads.
				is_read1_reversed = 0;
				is_read2_reversed = 0;
				strcpy(read_text_1, read_text_1_raw);
				strcpy(read_text_2, read_text_2_raw);
				if(is_null_qual){
					qual_text_1[0]=0;
					qual_text_2[0]=0;
				}
				free(read_text_2_raw);
				free(read_text_1_raw);
			}
		}

	}

	free(read_text_1);
	free(read_text_2);
	free(qual_text_1);
	free(qual_text_2);
	free(read_name_1);
	free(read_name_2);
	free(output_line_buffer);
	destroy_output_context(global_context, &out_context);	

	return 0;
}



int write_chunk_results_144(global_context_t * global_context)
{
	unsigned int read_number, sqr_read_number = 0, sqr_interval;
	gene_input_t * ginp1, * ginp2=NULL;
	char * additional_information = malloc(1800);
	short current_display_offset = 0, current_display_tailgate = 0;

	unsigned int out_poses[CIGAR_PERFECT_SECTIONS+1], xk1;
	char * out_cigars[CIGAR_PERFECT_SECTIONS+1], *out_mate_cigars[CIGAR_PERFECT_SECTIONS+1];
	char out_strands[CIGAR_PERFECT_SECTIONS+1];
	short out_read_lens[CIGAR_PERFECT_SECTIONS+1];

	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS+1; xk1++) out_cigars[xk1] = malloc(100);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS+1; xk1++) out_mate_cigars[xk1] = malloc(100);

	//if(global_context -> config.space_type == GENE_SPACE_COLOR && !global_context -> config.convert_color_to_base)
	//	current_display_offset = 1;


	ginp1 = &global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		ginp2 = &global_context->input_reads.second_read_file;
	
	sqr_interval = global_context -> processed_reads_in_chunk/10; 

	for(read_number = 0; read_number < global_context -> processed_reads_in_chunk ; read_number++)
	{
		char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
		char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
		char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
		int read_len_1, read_len_2=0;
		int best_read_id = 0;
		int best_read_id_HI = 0; 
		int total_best_reads = 0;
		int read1_has_reversed = 0;
		int read2_has_reversed = 0;
	
		int is_second_read;
		int applied_reverse_space;

		if(sqr_read_number > sqr_interval)
		{
			show_progress(global_context, NULL, read_number, STEP_WRITE_CHUNK_RESULTS);
			sqr_read_number = 0;
		}

		sqr_read_number ++;
		fetch_next_read_pair(global_context, NULL, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 0);
		//printf("ORAW1=%s\nORAW2=%s\n\n", read_text_1,read_text_2);

		applied_reverse_space = global_context->config.space_type;
		if(global_context -> config.convert_color_to_base)
		{
			colorread2base(read_text_1, read_len_1+1);
			if(global_context->input_reads.is_paired_end_reads)
				colorread2base(read_text_2, read_len_2+1);
			applied_reverse_space = GENE_SPACE_BASE;

		}
		//printf("BAS1=%s\nBAS2=%s\n\n", read_text_1,read_text_2);

		//alignment_result_t * prime_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0, 0);
		int read_1_repeats = 0, read_1_reported=0;
		int read_2_repeats = 0, read_2_reported=0;
		//int is_PE_OK = -1;
		for(total_best_reads=0; total_best_reads<global_context -> config.multi_best_reads; total_best_reads++)
		{
			alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context  , read_number, 0, total_best_reads);
			alignment_result_t * current_result_2 = NULL;

			//if(total_best_reads > 0 && current_result->selected_votes < 1) break;
			//if(total_best_reads > 0 && global_context -> input_reads.is_paired_end_reads && !is_result_in_PE(prime_result)) break;
			int read_1_result_available = 1;
			int read_2_result_available = 1;

			if(global_context -> config.multi_best_reads>1)
			{
				char * chro1_chro = NULL;
				unsigned int chro_1_pos = 0;
				if(locate_gene_position_max(current_result -> selected_position, &global_context -> chromosome_table, &chro1_chro, &chro_1_pos, read_len_1))
				{
					read_1_result_available = 0;
				}
			}

			if(read_1_result_available)
				if((current_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
				{
					read_1_result_available = 0;
				//	SUBREADprintf("Disabled read on time UNIQ:%s\n",current_name);
				}

			if(global_context -> input_reads.is_paired_end_reads)
			{
				current_result_2 = _global_retrieve_alignment_ptr(global_context  , read_number, 1, total_best_reads);

				if(read_2_result_available)
					if((current_result_2->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
					{
						read_2_result_available = 0;
					}


				if(global_context -> config.multi_best_reads>1)
				{
					char * chro1_chro = NULL;
					unsigned int chro_1_pos = 0;
					if(locate_gene_position_max(current_result_2 -> selected_position, &global_context -> chromosome_table, &chro1_chro, &chro_1_pos, read_len_2))
					{
						read_2_result_available = 0;
					}
				}

				if(current_result_2 -> result_flags & CORE_IS_FULLY_EXPLAINED)if(read_2_result_available) read_2_repeats ++;
				if(current_result -> result_flags & CORE_IS_FULLY_EXPLAINED)if(read_1_result_available) read_1_repeats ++;
			}else
			{
				if(current_result -> result_flags & CORE_IS_FULLY_EXPLAINED)
					if(read_1_result_available)read_1_repeats ++;
			}
		}

		for(best_read_id=0; best_read_id<global_context -> config.multi_best_reads; best_read_id++)
		{
			char mate_cigar_decompress[100];
			char current_cigar_decompress[100];
			for(is_second_read = 0; is_second_read < 1+ global_context -> input_reads.is_paired_end_reads; is_second_read++)
			{
				alignment_result_t *current_result, *mate_result = NULL;
				current_result = _global_retrieve_alignment_ptr(global_context  , read_number, is_second_read, best_read_id);
				
				if(global_context -> input_reads.is_paired_end_reads)
					mate_result    = _global_retrieve_alignment_ptr(global_context  , read_number,!is_second_read, best_read_id);

				//if(best_read_id>0 && current_result->selected_votes < 1 && (is_second_read == 0 || !global_context -> input_reads.is_paired_end_reads)) break;
				//if(best_read_id>0 && ( global_context -> input_reads.is_paired_end_reads&& !is_result_in_PE(prime_result))) break;
				/*
				if(global_context -> input_reads.is_paired_end_reads)
				{
					if(best_read_id)
						if((!(current_result -> result_flags &CORE_IS_FULLY_EXPLAINED)) && !(CORE_IS_FULLY_EXPLAINED &  mate_result -> result_flags))
							break;
				}else
				{
					if(best_read_id)
						if(!(current_result -> result_flags &CORE_IS_FULLY_EXPLAINED))
							break;
				}*/

				char * current_name      = is_second_read ? read_name_2 : read_name_1;
				char * current_read_text = is_second_read ? read_text_2 : read_text_1;
				char * current_qual_text = is_second_read ? qual_text_2 : qual_text_1;
				int * current_has_reversed = is_second_read ? (&read2_has_reversed):(&read1_has_reversed);
				int current_read_len     = is_second_read ? read_len_2  : read_len_1;
				int mate_read_len	= is_second_read ? read_len_1  : read_len_2;
				unsigned int mate_linear_pos=0, current_linear_pos=0;
				char mate_strand = 0, current_strand = 0;
				//int current_read_repeats = is_second_read ? read_2_repeats:read_1_repeats;
				int * current_read_reported = is_second_read ? (&read_2_reported):(&read_1_reported);

				char * current_chro_name=NULL, * mate_chro_name = NULL;
				unsigned int current_chro_offset=0, mate_chro_offset=0;
				char * current_CIGAR;
				int second_char = -1;

				additional_information[0]=0;

				int mask = 0;
				int is_mate_ok = 0;
				int is_current_ok = (current_result -> result_flags & CORE_IS_FULLY_EXPLAINED)?1:0;
				int current_repeated_times = 0;
				float current_final_quality = current_result -> final_quality;
				int current_soft_clipping_movement  =0, mate_soft_clipping_movement = 0;

				current_linear_pos = current_result -> selected_position;
				current_strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(global_context -> input_reads.is_paired_end_reads)
				{
					remove_backslash(current_name);

					mask |= SAM_FLAG_PAIRED_TASK;
					is_mate_ok = (mate_result -> result_flags & CORE_IS_FULLY_EXPLAINED)?1:0;

					if((mate_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
						is_mate_ok = 0;

					int mate_repeated_times = 0;

					mate_strand = (mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
					mate_linear_pos = mate_result -> selected_position;

					if(global_context -> config.do_big_margin_filtering_for_reads)
						mate_repeated_times = is_ambiguous_voting(global_context, read_number, !is_second_read, mate_result->selected_votes, mate_result->confident_coverage_start, mate_result->confident_coverage_end, mate_read_len, mate_strand);

					if( global_context -> config.do_big_margin_filtering_for_reads && mate_repeated_times > 1)
						is_mate_ok = 0;


					if(global_context->config.report_no_unpaired_reads && !is_result_in_PE(current_result))
					{
						is_mate_ok = mate_result->selected_votes > current_result->selected_votes?is_mate_ok:0;
						is_current_ok = mate_result->selected_votes < current_result->selected_votes?is_current_ok:0;
					}


					if(is_second_read)
						mask |= SAM_FLAG_SECOND_READ_IN_PAIR;
					else
						mask |= SAM_FLAG_FIRST_READ_IN_PAIR;



					if(is_mate_ok)
					{

						int is_jumped = 0;
						char * mate_CIGAR;
						if( mate_result -> cigar_string[0] == -1)
						{
							is_jumped = 1;
							bincigar2cigar(mate_cigar_decompress, 100, mate_result -> cigar_string + 1, CORE_MAX_CIGAR_LEN - 1, mate_read_len);

							mate_CIGAR = mate_cigar_decompress;
						}
						else
						{
							bincigar2cigar(mate_cigar_decompress, 100, mate_result -> cigar_string, CORE_MAX_CIGAR_LEN, mate_read_len);
							mate_CIGAR = mate_cigar_decompress;

						}

						if(global_context -> config.do_fusion_detection)
						{
							chimeric_cigar_parts(global_context, mate_linear_pos, is_jumped ^ mate_strand, is_jumped , mate_cigar_decompress , out_poses, out_mate_cigars, out_strands, mate_read_len, out_read_lens);

							mate_linear_pos = out_poses[0];
							mate_strand = out_strands[0]=='-';
							mate_CIGAR = out_mate_cigars[0];
						}

						if(locate_gene_position_max(mate_linear_pos, &global_context -> chromosome_table, & mate_chro_name, & mate_chro_offset, mate_read_len))
						{
							is_mate_ok = 0;
							//if(!is_second_read)
							//	read_2_repeats--;
								
						}
						mate_soft_clipping_movement = get_soft_clipping_length(mate_CIGAR);
						mate_chro_offset += mate_soft_clipping_movement;

					}
					if(is_mate_ok)
					{
						if(mate_strand + (!is_second_read) == 1) mask |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						mate_chro_offset++;
		
					}

				}
				if(!is_mate_ok)
				{
					mate_chro_name = "*";
					mate_chro_offset = 0;
				}
		

				if(global_context -> config.do_big_margin_filtering_for_reads)
					current_repeated_times = is_ambiguous_voting(global_context, read_number, is_second_read, current_result->selected_votes, current_result->confident_coverage_start, current_result->confident_coverage_end, current_read_len, current_strand);
		


				if(is_current_ok)
					if((current_result->result_flags & CORE_IS_BREAKEVEN) && !global_context -> config.report_multi_mapping_reads)
					{
						is_current_ok = 0;
					//	SUBREADprintf("Disabled read on time UNIQ:%s\n",current_name);
					}

				if(is_current_ok)
				{
					int is_first_section_jumped = 0;
					if( current_result -> cigar_string[0] == -1)
					{
						bincigar2cigar(current_cigar_decompress, 100, current_result -> cigar_string + 1, CORE_MAX_CIGAR_LEN, current_read_len);
						//current_linear_pos = reverse_cigar(current_linear_pos , current_cigar_decompress, current_cigar_decompress_new);
						current_CIGAR  = current_cigar_decompress;
						is_first_section_jumped = 1;
					}
					else
					{
						bincigar2cigar(current_cigar_decompress, 100, current_result -> cigar_string, CORE_MAX_CIGAR_LEN, current_read_len);
						current_CIGAR  = current_cigar_decompress;
					}


					if(global_context -> config.do_fusion_detection)
					{

						int chimeric_sections = chimeric_cigar_parts(global_context, current_linear_pos, is_first_section_jumped ^ current_strand, is_first_section_jumped , current_CIGAR , out_poses, out_cigars, out_strands, current_read_len, out_read_lens);

						//sprintf(additional_information + strlen(additional_information), "\tXX:Z:%s", current_cigar_decompress);

						for(xk1=1; xk1<chimeric_sections; xk1++)
						{
							unsigned int chimeric_pos;
							char * chimaric_chr;

							if(0==locate_gene_position_max(out_poses[xk1],& global_context -> chromosome_table, & chimaric_chr, & chimeric_pos, 0+out_read_lens[xk1]))
							{
								int soft_clipping_movement = 0;
								soft_clipping_movement = get_soft_clipping_length( out_cigars[xk1]);
								char strand_xor = (out_strands[xk1] == '-')^ is_second_read;
								sprintf(additional_information + strlen(additional_information), "\tCG:Z:%s\tCP:i:%u\tCT:Z:%c\tCC:Z:%s", out_cigars[xk1] , chimeric_pos + soft_clipping_movement + 1, strand_xor?'-':'+' , chimaric_chr );
							}
						}
							

						current_linear_pos = out_poses[0];
						current_strand = out_strands[0]=='-';
						current_CIGAR = out_cigars[0];
					}

					if(locate_gene_position_max(current_linear_pos,& global_context -> chromosome_table, & current_chro_name, & current_chro_offset, current_read_len))
					{
						is_current_ok = 0;
					//	if (!is_second_read) 
					//		read_1_repeats--;
					//	current_read_repeats--;
					}
				}

				if(is_current_ok)
				{
					if(current_strand + is_second_read == 1)
						mask |= SAM_FLAG_REVERSE_STRAND_MATCHED;
					if(*current_read_reported){
						mask |= SAM_FLAG_SECONDARY_MAPPING;
					}
					(*current_read_reported)=1;
					//if(1639 == read_number)
					//	printf("R0=%d ; NEG=%d; SEC=%d\n",(*current_has_reversed), (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0, is_second_read);
					int current_need_reverse = current_strand;

					//if(current_result -> cigar_string[0]==-1)
					//	current_need_reverse = ! current_need_reverse;

					if(current_need_reverse + (*current_has_reversed) == 1)
					{
						reverse_read(current_read_text, current_read_len + global_context->config.convert_color_to_base, applied_reverse_space);
						reverse_quality(current_qual_text , current_read_len);
						(*current_has_reversed)=!(*current_has_reversed);
					//	if(1639 == read_number)
					//		printf("RR=%d\n\n",(*current_has_reversed));
					}
					current_chro_offset++;

					current_soft_clipping_movement = get_soft_clipping_length(current_CIGAR);
					current_chro_offset += current_soft_clipping_movement;
				}
				else
				{
					mask |= SAM_FLAG_UNMAPPED;
					int this_should_nagetive = is_second_read;

					if(global_context -> input_reads.is_paired_end_reads && global_context -> config.report_unmapped_using_mate_pos&& is_mate_ok)
					{

						// DO NOT SHOW CORRDINATE IF IT IS UNMAPPED.
						//
						//current_chro_name = mate_chro_name;
						//current_chro_offset = mate_chro_offset;
						current_chro_name = "*";
						current_chro_offset = 0;

						/////////////////////////////////////////

						this_should_nagetive = mate_strand;
						current_strand = mate_strand;
						if(this_should_nagetive + is_second_read ==1)
							mask |= SAM_FLAG_REVERSE_STRAND_MATCHED;
						else
							mask &= ~SAM_FLAG_REVERSE_STRAND_MATCHED;
					}
					else
					{
						current_chro_name = "*";
						current_chro_offset = 0;
					}

					

					if(this_should_nagetive + (*current_has_reversed) == 1)
					{
						reverse_read(current_read_text, current_read_len + global_context->config.convert_color_to_base, applied_reverse_space);
						reverse_quality(current_qual_text , current_read_len);
						(*current_has_reversed)=!(*current_has_reversed);
					}

					current_CIGAR = "*";
					current_final_quality=0;
				}

				if(global_context -> config.space_type == GENE_SPACE_COLOR)
				{
					if( is_second_read  +  current_strand == 1 )
					{
						//if(is_current_ok)
						//	current_chro_offset ++;
						current_display_offset = 0;
						current_display_tailgate = 1;
					}
					else
					{
						if(!global_context -> config.convert_color_to_base)
						{
							// the first base was a fake prime base; the second base is the first meaningful base.
							second_char = current_read_text[1];
							current_read_text[1] = color2char(second_char, current_read_text[0]);
						}
						current_display_offset = 1;
						current_display_tailgate = 0;
					}
				}

				if(global_context -> input_reads.is_paired_end_reads && !is_mate_ok)
				{
					mask |= SAM_FLAG_MATE_UNMATCHED;
					mate_strand = current_strand;

					if(is_current_ok){

						int mate_should_nagetive = mate_strand;
						if(mate_should_nagetive + (!is_second_read) ==1)
							mask |= SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
						else
							mask &= ~SAM_FLAG_MATE_REVERSE_STRAND_MATCHED;
					}


				}
				long long int mate_distance = 0;

				// DO NOT SHOW CORRDINATE IF IT IS UNMAPPED.
				//
					if( 0 && is_current_ok && global_context -> config.report_unmapped_using_mate_pos && global_context -> input_reads.is_paired_end_reads &&!is_mate_ok)
					{
						mate_distance = 0;

						
						mate_chro_name = current_chro_name;
						mate_chro_offset = current_chro_offset;
					}
				//////////////////////////////////////////////


				if(is_current_ok && global_context -> input_reads.is_paired_end_reads && is_mate_ok)
				{
					mate_distance = mate_chro_offset - mate_soft_clipping_movement;
					mate_distance -= current_chro_offset - current_soft_clipping_movement;
					mate_distance = abs(mate_distance);
					if(current_chro_offset >mate_chro_offset)
						mate_distance += current_read_len; 
					else
						mate_distance += mate_read_len; 

					if(current_chro_offset - current_soft_clipping_movement > mate_chro_offset - mate_soft_clipping_movement) mate_distance = -mate_distance;

					if(mate_distance>0)
					{
						mate_distance = max(mate_distance, current_read_len);
						mate_distance = max(mate_distance, mate_read_len);
					}
					else
					{
						mate_distance = min(mate_distance, -current_read_len);
						mate_distance = min(mate_distance, -mate_read_len);
					}
				}

				if((best_read_id == 0 ) && current_qual_text[0] && global_context -> config.phred_score_format == FASTQ_PHRED64)
					fastq_64_to_33(current_qual_text);
				if(!current_qual_text[0])
				{
					int xi2;
					for(xi2=current_display_offset;current_read_text[xi2 + current_display_tailgate];xi2++) current_qual_text[xi2 - current_display_offset] = 'I';
					current_qual_text[xi2 - current_display_offset]=0;
				}

				if(current_repeated_times)
					current_final_quality /= current_repeated_times;
				if(global_context->config.downscale_mapping_quality)
					current_final_quality=current_final_quality/5;

				if(mate_chro_name[0]!='*' && mate_chro_name == current_chro_name && global_context -> input_reads.is_paired_end_reads)
				{
					mate_chro_name="=";

					//if(2190 == read_number)SUBREADprintf("MATE_DIST=%lld\n", mate_distance);

					if(is_mate_ok && is_current_ok && abs(mate_distance) >= global_context->config. minimum_pair_distance && abs(mate_distance) <= global_context-> config.maximum_pair_distance  && current_strand == mate_strand)
					{
					//	if(2190 == read_number)SUBREADprintf("MATE_DIST 2 =%lld , CUR_OFF=%u, MAT_OFF=%u, CUR_STRAND=%d\n", mate_distance, current_chro_offset, mate_chro_offset, current_strand);
						int is_PEM = 0;
						if(global_context -> config.do_fusion_detection)
						{
							is_PEM = 1;
						}
						else 
						{
							if(global_context -> config.is_first_read_reversed && !(global_context -> config.is_second_read_reversed))
							{
								if(current_strand == 0)
								{
									if((is_second_read + (mate_chro_offset > current_chro_offset) == 1) || mate_chro_offset == current_chro_offset)
										is_PEM = 1;
								}
							}
							else
							{
								if(current_strand)
								{
									if((is_second_read + (mate_chro_offset < current_chro_offset) == 1) || mate_chro_offset == current_chro_offset) is_PEM = 1;
								}else{
									if((is_second_read + (mate_chro_offset > current_chro_offset) == 1) || mate_chro_offset == current_chro_offset) is_PEM = 1;

								}
							}
						}

						//if(2190 == read_number)SUBREADprintf("MATE_DIST 3 =%lld ; PEM=%d\n", mate_distance, is_PEM);
						if(is_PEM){
							mask |= SAM_FLAG_MATCHED_IN_PAIR;
							if(is_second_read)
								global_context->all_correct_PE_reads  ++;
						}
					}


				}
				if(mate_chro_name[0]!='=') 
					mate_distance = 0;

				int tailgate_0 = -1;
				if(current_display_tailgate)
				{
					tailgate_0 = current_read_text[strlen(current_read_text) -1];
					current_read_text[strlen(current_read_text) - 1] = 0;
				}

				//if(read_number==2461)
				//printf("CURCIGAR=%s ; FINAL_CIGAR=%s ; OK=%d\n", current_cigar_decompress, current_CIGAR, is_current_ok);

				if( is_mate_ok ||is_current_ok || (((global_context -> config.multi_best_reads <2) || (best_read_id==0 && read_1_repeats == 0 && read_2_repeats == 0)) && ! global_context->config.ignore_unmapped_reads))
				{

					char hi_tag_out[18];
					hi_tag_out[0]=0;

					if(max(read_2_repeats,read_1_repeats) > 1)
						sprintf(hi_tag_out,"\tHI:i:%d", best_read_id_HI);


					int seq_len = strlen(additional_information);
					seq_len += sprintf(additional_information+seq_len, "\tSH:i:%d\tSM:i:%d\tNH:i:%d%s", (int)((current_result -> Score_H >> 17) & 0xfff), current_result -> final_mismatched_bases, max(read_2_repeats,read_1_repeats), hi_tag_out);

					if( is_current_ok && global_context -> config.is_rna_seq_reads && !(current_result -> result_flags & CORE_NOTFOUND_DONORS))
					{
						seq_len += sprintf(additional_information+seq_len, "\tXS:A:%c", (current_result -> result_flags & CORE_IS_GT_AG_DONORS)?'+':'-');
					}

					if(global_context -> config.SAM_extra_columns)
					{
						if(!is_current_ok){
							current_cigar_decompress[0]='*';
							current_cigar_decompress[1]=0;
						}
						seq_len += sprintf(additional_information+seq_len, "\tSG:Z:%s\tSB:i:%d\tSC:i:%d\tSD:i:%d\tSN:i:%u\tSP:Z:%s",current_cigar_decompress, current_result -> used_subreads_in_vote, current_result -> selected_votes, current_result -> noninformative_subreads_in_vote, read_number, (current_result -> result_flags & CORE_IS_GAPPED_READ)?"GAPPED":"NOEVENT"); 
					}

					if(global_context->config.read_group_id[0])
						seq_len += sprintf(additional_information+seq_len, "\tRG:Z:%s", global_context->config.read_group_id);

					if(global_context -> config.is_BAM_output)
						SamBam_writer_add_read(global_context -> output_bam_writer, current_name, mask, current_chro_name, current_chro_offset, (int)current_final_quality, current_CIGAR, mate_chro_name, mate_chro_offset, mate_distance, current_read_len, current_read_text + current_display_offset, current_qual_text, additional_information+1);
					else
						sambamout_fprintf(global_context -> output_sam_fp , "%s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%s\t%s%s\n", current_name, mask, current_chro_name, current_chro_offset, (int)current_final_quality, current_CIGAR, mate_chro_name, mate_chro_offset, mate_distance, current_read_text + current_display_offset, current_qual_text, additional_information);
				}

				if(current_display_tailgate)
				{
					current_read_text[strlen(current_read_text)] = tailgate_0;
				}

				if(second_char > 0)
					current_read_text[1] = second_char;

				if(global_context->input_reads.is_paired_end_reads)
				{
					if(is_second_read)
						if(is_current_ok || is_mate_ok)
						{
							global_context -> all_mapped_reads++;
							best_read_id_HI++;
						}
				}
				else
				{
					if(is_current_ok)
					{
						best_read_id_HI++;
						global_context -> all_mapped_reads++;
					}
				}
			} 
		}
	}

	free(additional_information);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS; xk1++) free(out_cigars[xk1]);
	for(xk1 = 0; xk1 < CIGAR_PERFECT_SECTIONS; xk1++) free(out_mate_cigars[xk1]);
	return 0;
}
void init_chunk_scanning_parameters(global_context_t * global_context, thread_context_t * thread_context, gene_input_t ** ginp1, gene_input_t ** ginp2, unsigned int * read_block_start, unsigned int * reads_to_be_done)
{
	*ginp2 = NULL;
	*ginp1 = thread_context?thread_context->ginp1: & global_context->input_reads.first_read_file;
	if(global_context->input_reads.is_paired_end_reads)
		*ginp2 = thread_context?thread_context->ginp2:& global_context->input_reads.second_read_file;

	*read_block_start = thread_context?thread_context->read_block_start:0;
	*reads_to_be_done = thread_context?thread_context->reads_to_be_done:global_context -> config.reads_per_chunk;
}

gene_value_index_t * find_current_value_index(global_context_t * global_context, unsigned int pos, int len)
{
	int block_no;
	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		if(pos>=index_begin && pos + len<index_end)
			return & global_context->all_value_indexes [0];
		else return NULL;
	}
	else
		for(block_no=0; block_no<global_context->index_block_number; block_no++)
		{
			unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
			unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
			if((block_no == 0 && pos >=  index_begin && pos < index_end - 1000000) ||
			   (block_no>0 && block_no<global_context->index_block_number-1 && pos >= index_begin+ 1000000 && pos < index_end - 1000000) ||
			   (block_no == global_context->index_block_number-1 && pos >= index_begin + 1000000 && pos < index_end ))
			{
				return & global_context -> all_value_indexes [block_no];
			}
		}
	return NULL;
}
//this function selects the correct all_value_indexes from global_context and put it to global_context or thread_context if thread_context is not NULL.
int locate_current_value_index(global_context_t * global_context, thread_context_t * thread_context, alignment_result_t * result, int rlen)
{
	int block_no;

	if(global_context->index_block_number<2)
	{
		unsigned index_begin = global_context -> all_value_indexes [0] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [0] . start_base_offset + global_context -> all_value_indexes [0] . length;

		//SUBREADprintf("RESET2 : %u should <= %u\n",  result->selected_position + rlen, index_end);

		if(result->selected_position>=index_begin && result->selected_position + rlen<=index_end)
		{
			if(thread_context)thread_context->current_value_index = & global_context->all_value_indexes [0];
			else global_context->current_value_index =& global_context->all_value_indexes [0];
			return 0;
		}
		else return 1;
	}
	for(block_no=0; block_no<global_context->index_block_number; block_no++)
	{
		unsigned index_begin = global_context -> all_value_indexes [block_no] . start_base_offset; 
		unsigned index_end = global_context -> all_value_indexes [block_no] . start_base_offset + global_context -> all_value_indexes [block_no] . length;
		if((block_no == 0 && result->selected_position >=  index_begin && result->selected_position < index_end - 1000000) ||
		   (block_no>0 && block_no<global_context->index_block_number-1 && result->selected_position >= index_begin+ 1000000 && result->selected_position < index_end - 1000000) ||
		   (block_no == global_context->index_block_number-1 && result->selected_position >= index_begin + 1000000 && result->selected_position < index_end ))
		{
			if(thread_context)thread_context->current_value_index =& global_context -> all_value_indexes [block_no];
			else global_context->current_value_index =& global_context -> all_value_indexes [block_no];
			return 0;
		}
	}
	return 1;
}




int do_iteration_one(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	int need_junction_step = global_context -> config.is_rna_seq_reads || global_context -> config.do_fusion_detection;
	int sqr_interval, sqr_read_number = 0;

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 1);
		// if no more reads
		if(ret)
			break;

		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int best_read_id, is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id); 

	//	if(strcmp("a4", read_name_1) == 0)
	//		printf("IDR=%u   VOT=%d  PAIR#=%u\n", current_result->selected_position, current_result->selected_votes, current_read_number);



				if(current_result -> selected_votes<1) break;
				if(!global_context->config.report_multi_mapping_reads)if(current_result -> result_flags & CORE_IS_BREAKEVEN) continue;

				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				char * current_read_name = is_second_read?read_name_2:read_name_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;

				if(current_result->selected_votes < global_context->config.minimum_subread_for_first_read)
					continue;
				int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(is_negative_strand + is_reversed_already == 1)
				{
					reverse_read(current_read, current_rlen ,  global_context->config.space_type);
					reverse_quality(current_qual, current_rlen);
					is_reversed_already = !is_reversed_already;
				}

				if(locate_current_value_index(global_context, thread_context, current_result, current_rlen))
				{
				//	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR, "Read position excesses index boundary : %u (%s : %s). V=%d", current_result -> selected_position, current_read_name, is_second_read?"SECOND":"FIRST", current_result -> selected_votes);
					continue;
				}

				find_new_indels(global_context, thread_context, current_read_number, current_read_name, current_read, current_qual, current_rlen, is_second_read, best_read_id);
				if(need_junction_step)
					find_new_junctions(global_context, thread_context, current_read_number, current_read, current_qual, current_rlen, is_second_read, best_read_id);
			}
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number > sqr_interval)	
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_ONE);
				sqr_read_number = 0;
			}
		}

	}



	return 0;
}



int finish_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	return 0;
}
int do_iteration_three(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0, sqr_interval, sqr_read_number = 0;

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	print_in_box(80,0,0,"Prepare for long indel deleteion...");
	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2, 1);
		if(ret)
			break;

		int best_read_id = 0;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t * current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);
				alignment_result_t * mate_result   = _global_retrieve_alignment_ptr(global_context, current_read_number,!is_second_read, best_read_id); 

				if(best_read_id && (current_result->selected_votes <1)) break;

				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int mate_rlen = is_second_read?read_len_1:read_len_2;

				//if(!is_ambiguous_voting(global_context , current_read_number, is_second_read, current_result -> selected_votes, current_result -> confident_coverage_start, current_result -> confident_coverage_end, current_rlen, (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0))
				{
					// do local reassambly
					// a potential long-indel read has to have minimum supporting subreads, but not as many as total_subread - 1

					if(current_result->selected_votes >= global_context->config.minimum_subread_for_first_read)
					{
						int is_negative_strand = ((current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);

						if(is_negative_strand + is_reversed_already ==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number, current_read_name , current_read , current_qual , current_rlen, 0 , is_second_read, best_read_id, 0);

					}
					else if(global_context -> input_reads.is_paired_end_reads && is_result_in_PE(current_result) && current_result -> selected_votes >= global_context->config.minimum_subread_for_second_read)
					{
						int is_negative_strand = (current_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;
						if(is_negative_strand + is_reversed_already==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number, current_read_name , current_read , current_qual , current_rlen , 0, is_second_read, best_read_id, 0);
					}
					else if(global_context -> input_reads.is_paired_end_reads && mate_result -> selected_votes >= global_context->config.minimum_subread_for_first_read)
					{
						int is_negative_strand = ((mate_result -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0);
						if(is_negative_strand+is_reversed_already==1)
						{
							reverse_read(current_read, current_rlen,  global_context->config.space_type);
							reverse_quality(current_qual, current_rlen);
							is_reversed_already=!is_reversed_already;
						}

						build_local_reassembly(global_context , thread_context , current_read_number , current_read_name , current_read , current_qual , current_rlen, mate_rlen , is_second_read, best_read_id, 1);
					}
				}	
			}
		}

		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_THREE);
				sqr_read_number = 0;
			}
		}
	}

	return 0;

}


int do_iteration_two(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char read_text_1[MAX_READ_LENGTH+1], read_text_2[MAX_READ_LENGTH+1];
	char qual_text_1[MAX_READ_LENGTH+1], qual_text_2[MAX_READ_LENGTH+1];
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	int sqr_interval, sqr_read_number=0;

	//unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	//unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);
	sqr_interval = global_context -> processed_reads_in_chunk/10/ global_context -> config.all_threads;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;

		sqr_read_number++;
		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2,1);
		// if no more reads
		if(ret)
			break;

		int best_read_id=0;
		for (is_second_read = 0; is_second_read < 1 + global_context -> input_reads.is_paired_end_reads; is_second_read ++)
		{
			int is_reversed_already = 0;
			for(best_read_id = 0; best_read_id < global_context -> config.multi_best_reads; best_read_id++)
			{
				alignment_result_t *current_result = _global_retrieve_alignment_ptr(global_context, current_read_number, is_second_read, best_read_id);
				if(best_read_id > 0 && current_result -> selected_votes==0)break;
				if(!global_context->config.report_multi_mapping_reads)if(current_result -> result_flags & CORE_IS_BREAKEVEN) continue;

				char * current_read_name =  is_second_read?read_name_2 : read_name_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;


				if(current_result->selected_votes < global_context->config.minimum_subread_for_second_read)
				{
					//SUBREADprintf("RESET0 %d %d\n", current_result->selected_votes);
					current_result->selected_votes = 0;
					current_result -> final_mismatched_bases = 0;
					continue;
				}
				if(locate_current_value_index(global_context, thread_context, current_result, current_rlen))
				{
					current_result->selected_votes = 0;
					current_result -> final_mismatched_bases = 0;
					//SUBREADprintf("RESET1 Read position excesses index boundary.\n");
					continue;
				}

				int is_negative_strand = (current_result  -> result_flags & CORE_IS_NEGATIVE_STRAND)?1:0;

				if(is_negative_strand + is_reversed_already == 1)
				{
					reverse_read(current_read, current_rlen,  global_context->config.space_type);
					reverse_quality(current_qual, current_rlen);
					is_reversed_already = !is_reversed_already;
				}

				explain_read(global_context, thread_context, current_read_number, current_rlen, current_read_name, current_read, current_qual, is_second_read, best_read_id, is_negative_strand);
			}
		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number>sqr_interval)
			{
				show_progress(global_context, thread_context, current_read_number, STEP_ITERATION_TWO);
				sqr_read_number=0;
			}
		}
	}

	return 0;
}

int core_get_subread_quality(global_context_t * global_context, thread_context_t * thread_context, char * qual, int qual_len)
{
	int x1, ret=1;

	if(!qual)return 1;
	if(!qual[0])return 1;

	int offset = (global_context->config.phred_score_format == FASTQ_PHRED33)?33:64; 

	for(x1=0; (x1 < qual_len) && qual[x1]; x1++)
		ret +=  (qual[x1] - offset);

	return  ret;
}

int do_voting(global_context_t * global_context, thread_context_t * thread_context)
{
	unsigned int reads_to_be_done = 0, read_block_start = 0;
	int ret, xk1;
	gene_input_t * ginp1 = NULL , * ginp2 = NULL;
	unsigned int current_read_number;
	char * read_text_1, * read_text_2;
	char * qual_text_1, * qual_text_2;
	char read_name_1[MAX_READ_NAME_LEN+1], read_name_2[MAX_READ_NAME_LEN+1];
	int read_len_1, read_len_2=0;
	unsigned int processed_reads=0;
	int min_first_read_votes = global_context -> config.minimum_subread_for_first_read; 
	int voting_max_indel_length = min(16, global_context->config.max_indel_length);
	int sqr_interval=10000, sqr_read_number = 0;

	read_text_1 = malloc(MAX_READ_LENGTH+1);
	read_text_2 = malloc(MAX_READ_LENGTH+1);
	qual_text_1 = malloc(MAX_READ_LENGTH+1);
	qual_text_2 = malloc(MAX_READ_LENGTH+1);

	gene_vote_t * vote_1, * vote_2, * vote_fg;

	vote_1 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_2 = (gene_vote_t *) malloc(sizeof(gene_vote_t));
	vote_fg = (gene_vote_t *) malloc(sizeof(gene_vote_t));

	if(!(vote_1&&vote_2))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_FATAL,"Cannot allocate voting memory.");
		return 1;
	}

	init_chunk_scanning_parameters(global_context,thread_context, & ginp1, & ginp2, & read_block_start, & reads_to_be_done);

	unsigned int low_index_border = global_context -> current_value_index -> start_base_offset;
	unsigned int high_index_border = global_context -> current_value_index -> start_base_offset + global_context -> current_value_index -> length; 
	int has_second_read = 1 + global_context -> input_reads.is_paired_end_reads;
	//int need_junction_step = global_context -> config.is_rna_seq_reads || global_context -> config.do_fusion_detection;

	if(thread_context)
		thread_context -> current_value_index = global_context -> current_value_index;

	int GENE_SLIDING_STEP = global_context -> current_index -> index_gap;

	for(current_read_number = read_block_start; current_read_number < reads_to_be_done + read_block_start ; current_read_number++)
	{
		int is_second_read;
		int subread_no;
		int is_reversed, applied_subreads = 0, v1_all_subreads=0, v2_all_subreads=0;

		ret = fetch_next_read_pair(global_context, thread_context, ginp1, ginp2, &read_len_1, &read_len_2, read_name_1, read_name_2, read_text_1, read_text_2, qual_text_1, qual_text_2,1);
		if(ret)
			break;

		//printf("%s\t%d\n%s\t%d\n", read_name_1, thread_context -> thread_id, read_name_2,  thread_context -> thread_id);

		for(is_reversed = 0; is_reversed<2; is_reversed++)
		{
			//printf("MAP_REA = %s / %s\n", read_text_1, read_text_2);
			if(is_reversed==0 || !global_context->config.do_fusion_detection)
			{
				init_gene_vote(vote_1);
				if(global_context -> input_reads.is_paired_end_reads) init_gene_vote(vote_2);
			}


			for (is_second_read = 0; is_second_read < has_second_read; is_second_read ++)
			{
				gene_vote_t * current_vote = is_second_read?vote_2: vote_1;
				char * current_read =  is_second_read?read_text_2 : read_text_1;
				char * current_qual =  is_second_read?qual_text_2 : qual_text_1;
				int current_rlen = is_second_read?read_len_2:read_len_1;
				int subread_step;
				if(current_rlen< 16) continue;
				if(current_rlen<= EXON_LONG_READ_LENGTH)
				{
					subread_step = (((current_rlen - 15 - GENE_SLIDING_STEP)<<16))/(global_context -> config.total_subreads -1);
					if(subread_step<(GENE_SLIDING_STEP<<16))subread_step = GENE_SLIDING_STEP<<16;
				}else{
					subread_step = 6<<16;
					if(((current_rlen - 15-GENE_SLIDING_STEP)<<16) / subread_step > 62)
						subread_step = ((current_rlen - 15-GENE_SLIDING_STEP)<<16)/62;
				}

				int allsubreads_for_each_gap [GENE_SLIDING_STEP], noninformative_subreads_for_each_gap[GENE_SLIDING_STEP];

				applied_subreads = 1 + ((current_rlen - 15-GENE_SLIDING_STEP)<<16) / subread_step;
				if(is_second_read) v2_all_subreads = applied_subreads;
				else	 v1_all_subreads = applied_subreads;

				//SUBREADprintf("NSUBR=%d\tAPPLIED_SUBR=%d\tSTEP=%d\n", global_context -> config.total_subreads, applied_subreads, subread_step);

				unsigned int current_high_border = high_index_border -  current_rlen;

				if(global_context->config.is_rna_seq_reads && current_rlen > EXON_LONG_READ_LENGTH && global_context->config.all_threads<2)
					core_fragile_junction_voting(global_context, thread_context, current_read, current_qual, current_rlen, is_reversed, global_context->config.space_type, low_index_border, current_high_border, vote_fg);

				if(global_context->config.SAM_extra_columns)
				{
					for(xk1=0;xk1<GENE_SLIDING_STEP;xk1++)
					{
						allsubreads_for_each_gap[xk1]=0;
						noninformative_subreads_for_each_gap[xk1]=0;
					}
				}

				for(subread_no=0; subread_no < applied_subreads ; subread_no++)
				{
					for(xk1=0; xk1<GENE_SLIDING_STEP ; xk1++)
					{

						if(global_context->config.SAM_extra_columns)
						{
							current_vote -> noninformative_subreads = noninformative_subreads_for_each_gap[xk1];
						}

						int subread_offset = ((subread_step * subread_no) >> 16);
						if(GENE_SLIDING_STEP > 1)
							subread_offset -= subread_offset%(GENE_SLIDING_STEP) - xk1; 

						int subread_quality = 1;
						char * subread_string = current_read + subread_offset;

						gehash_key_t subread_integer = genekey2int(subread_string, global_context->config.space_type);


						if(global_context -> config.use_quality_score_break_ties)
						{
							char * quality_string_subr = current_qual + subread_offset;
							subread_quality = core_get_subread_quality(global_context, thread_context, quality_string_subr, 16);
						}



						//SUBREADprintf("%d ", subread_offset);
						if(global_context->config.is_methylation_reads)
							gehash_go_q_CtoT(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, 1, subread_quality, 0xffffff, voting_max_indel_length, subread_no, 1,  low_index_border, high_index_border - current_rlen);
						else
							gehash_go_q(global_context->current_index, subread_integer , subread_offset, current_rlen, is_reversed, current_vote, subread_quality, voting_max_indel_length, subread_no,  low_index_border, current_high_border);



						if(global_context->config.SAM_extra_columns)
						{
							noninformative_subreads_for_each_gap[xk1] = current_vote -> noninformative_subreads;
						}


					}
					//SUBREADprintf(",");
				}
	//			SUBREADprintf("\n");

				//puts("");

				if(global_context->config.SAM_extra_columns)
				{
					short max_noninformative_subreads = -1;

					for(xk1=0;xk1<GENE_SLIDING_STEP;xk1++)
						if(noninformative_subreads_for_each_gap[xk1] > max_noninformative_subreads)
						{
							max_noninformative_subreads = noninformative_subreads_for_each_gap[xk1];
						}

					current_vote -> noninformative_subreads = max_noninformative_subreads;
				}
			}



			if(is_reversed==1 || !global_context->config.do_fusion_detection)
			{
				//if(strcmp(read_name_1,"b1")==0)

				//if(current_read_number == 119) {
				//	SUBREADprintf("NOINF=%d\n", vote_1 -> noninformative_subreads );
		//		#warning =============== COMMENT THIS LINE!!!! ======================
		//			print_votes(vote_1, global_context -> config.index_prefix);
				//}
				//	print_votes(vote_2, global_context -> config.index_prefix);
				//}

				//finalise_vote(vote_1);
 
				/*
				if(407229 == current_read_number){
					fprintf(stderr, "TABLE_ITEMS=%llu\n", global_context->current_index->current_items);
					print_votes(vote_1, global_context -> config.index_prefix);
					//print_votes(vote_2, global_context -> config.index_prefix);
				}*/
				//if(global_context -> input_reads.is_paired_end_reads) finalise_vote(vote_2);

				if(global_context -> input_reads.is_paired_end_reads)
					process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, read_name_2, read_text_1, read_text_2, read_len_1, read_len_2, is_reversed, v1_all_subreads, v2_all_subreads);
				else{
					if(vote_1->max_vote >= min_first_read_votes)
						process_voting_junction(global_context, thread_context, current_read_number, vote_1, vote_2, read_name_1, NULL ,  read_text_1, NULL, read_len_1, 0, is_reversed, v1_all_subreads, 0);
					else if(_global_retrieve_alignment(global_context, current_read_number, 0,0).selected_votes < 1)
					{
						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->used_subreads_in_vote = max(_global_retrieve_alignment(global_context, current_read_number, 0,0).used_subreads_in_vote, applied_subreads);
						_global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->noninformative_subreads_in_vote = max(_global_retrieve_alignment(global_context, current_read_number, 0,0).noninformative_subreads_in_vote, vote_1 -> noninformative_subreads);
					}
				}


				//if(current_read_number == 0) print_votes(vote_1, global_context -> config.index_prefix);
				//if(current_read_number == 0) printf("V0=%d; %d\n", _global_retrieve_alignment_ptr(global_context, current_read_number, 0,0)->selected_votes, _global_retrieve_alignment_ptr(global_context, current_read_number, 1,0)->selected_votes);
			}
		
			if(is_reversed == 0)
			{
				reverse_read(read_text_1, read_len_1,  global_context->config.space_type);
				reverse_quality(qual_text_1, read_len_1);

				if(global_context -> input_reads.is_paired_end_reads)
				{
					reverse_read(read_text_2, read_len_2,  global_context->config.space_type);
					reverse_quality(qual_text_2, read_len_2);
				}
			}

		}
		
		if(!thread_context || thread_context->thread_id == 0)
		{
			if(sqr_read_number > sqr_interval)
			{
				show_progress(global_context, thread_context, processed_reads, STEP_VOTING);
				sqr_read_number = 0;
				unsigned long long total_file_size = global_context -> input_reads.first_read_file_size;
				unsigned long long guessed_all_reads = total_file_size / global_context -> input_reads . avg_read_length;// / (1+global_context -> config.is_SAM_file_input);
				sqr_interval = guessed_all_reads / global_context -> config.all_threads/10;
			}
		}




		sqr_read_number++;
		processed_reads++;

	}

	if(thread_context)
		thread_context -> processed_reads_in_chunk = processed_reads;
	else
		global_context -> processed_reads_in_chunk = processed_reads;

	free(vote_1);
	free(vote_2);
	free(vote_fg);
	free(read_text_1);
	free(read_text_2);
	free(qual_text_1);
	free(qual_text_2);

	return 0;
}

void * run_in_thread(void * pthread_param)
{
	void ** parameters = (void **)pthread_param;
	global_context_t * global_context = (global_context_t * ) parameters[0];
	thread_context_t * thread_context = (thread_context_t *) parameters[1];
	int task = *((int *)(parameters[2]));
	subread_lock_t * thread_init_lock = (subread_lock_t * ) parameters[3];
	int * ret_value_pointer = (int *)parameters[4];

	if(thread_init_lock)
		subread_lock_release(thread_init_lock);

	switch(task)
	{
		case STEP_VOTING:
			*ret_value_pointer = do_voting(global_context, thread_context);
		break;
		case STEP_ITERATION_ONE:
			*ret_value_pointer = do_iteration_one(global_context, thread_context);
		break;
		case STEP_ITERATION_TWO:
			*ret_value_pointer = do_iteration_two(global_context, thread_context);
		break;
	
	}

	//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_DETAILS, "finished running %d", task);

	return NULL;
}

int run_maybe_threads(global_context_t *global_context, int task)
{
	void * thr_parameters [5];
	int ret_value =0;

	if(task==STEP_VOTING)
		print_in_box(80,0,0, "Map %s...", global_context->input_reads.is_paired_end_reads?"fragments":"reads");
	else if(task == STEP_ITERATION_ONE)
		print_in_box(80,0,0, "Detect indels%s...", global_context->config.is_rna_seq_reads?" and junctions":"");
	else if(task == STEP_ITERATION_TWO)
		print_in_box(80,0,0, "Realign %s...", global_context->input_reads.is_paired_end_reads?"fragments":"reads");

	if(global_context->config.all_threads<2)
	{
		thr_parameters[0] = global_context;
		thr_parameters[1] = NULL;
		thr_parameters[2] = &task;
		thr_parameters[3] = NULL;
		thr_parameters[4] = &ret_value;

		run_in_thread(thr_parameters);
	}
	else
	{
		int current_thread_no ;
		thread_context_t thread_contexts[64];
		int ret_values[64];

		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++)
		{
			thread_contexts[current_thread_no].thread_id = current_thread_no;
			init_indel_thread_contexts(global_context, thread_contexts+current_thread_no, task);
			if(global_context->config.is_rna_seq_reads || global_context->config.do_fusion_detection)
				init_junction_thread_contexts(global_context, thread_contexts+current_thread_no, task);

			relocate_geinputs(global_context, thread_contexts+current_thread_no);

			subread_lock_occupy(&global_context -> thread_initial_lock);
			thr_parameters[0] = global_context;
			thr_parameters[1] = thread_contexts+current_thread_no;
			thr_parameters[2] = &task;
			thr_parameters[3] = (void *)&global_context -> thread_initial_lock;
			thr_parameters[4] = ret_values + current_thread_no;

			pthread_create(&thread_contexts[current_thread_no].thread, NULL,  run_in_thread, &thr_parameters);
		}

		if(task == STEP_VOTING)
		{
			global_context -> processed_reads_in_chunk=0;
		}
		for(current_thread_no = 0 ; current_thread_no < global_context->config.all_threads ; current_thread_no ++)
		{
			pthread_join(thread_contexts[current_thread_no].thread, NULL);
			
			geinput_close(thread_contexts[current_thread_no].ginp1);
			free(thread_contexts[current_thread_no].ginp1);

			if(global_context->input_reads.is_paired_end_reads)
			{
				geinput_close(thread_contexts[current_thread_no].ginp2);
				free(thread_contexts[current_thread_no].ginp2);
			}

			ret_value += *(ret_values + current_thread_no);
			if(ret_value)break;

			if(task == STEP_VOTING)
			{
				//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_DEBUG, "The %d-th thread processed %u reads.", current_thread_no , thread_contexts[current_thread_no].processed_reads_in_chunk);
				global_context -> processed_reads_in_chunk += thread_contexts[current_thread_no].processed_reads_in_chunk;
			}

			finalise_indel_thread(global_context, thread_contexts+current_thread_no, task);
			finalise_junction_thread(global_context, thread_contexts+current_thread_no, task);
		}
	}

	if(CORE_SOFT_BR_CHAR == '\r')
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "");
	return ret_value;
}

void clean_context_after_chunk(global_context_t * context)
{
	memset(context -> chunk_alignment_records , 0 , sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
	memset(context -> big_margin_record  , 0 , sizeof(*context -> big_margin_record) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.big_margin_record_size);
	if(context ->chunk_subjunc_records)
		memset(context ->chunk_subjunc_records , 0 , sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
}

#define SKIP_CORE_NOEMPTY(fp_loc, buf_loc)	{ while(1){char *ret_loc = fgets(buf_loc, 3000, fp_loc); if(buf_loc[0]!='\n' || !ret_loc) break; } }

unsigned int split_read_files(global_context_t * global_context)
{
	unsigned int chunk_reads = global_context->config.reads_per_chunk;
	unsigned int processed_reads = 0;
	unsigned long long * read_position_1;
	unsigned long long * read_position_2 = NULL;
	char * read_line_buf = malloc(3002);
	char * read_line_buf2 = malloc(3002);
	read_position_1 = (unsigned long long*)malloc(global_context->config.reads_per_chunk * sizeof(long long));
	if(global_context->input_reads.is_paired_end_reads)
		read_position_2 = (unsigned long long*)malloc(global_context->config.reads_per_chunk * sizeof(long long));

	print_in_box(80,0,0, "Scan read files for multi-threaded alignment...");

	if(global_context->config.is_SAM_file_input)
	{
		unsigned long long fhead_pos1;
		unsigned long long fhead_pos2=0;

		while(1)
		{
			char * tok_tmp = NULL, * flag,  *flag2;
			if(processed_reads >= chunk_reads || feof(global_context->input_reads.first_read_file.input_fp))
				break;

			fhead_pos1 = ftello(global_context->input_reads.first_read_file.input_fp); 
			if(global_context->input_reads.is_paired_end_reads)
				fhead_pos2 = ftello(global_context->input_reads.second_read_file.input_fp); 

			fgets(read_line_buf, 3000, global_context->input_reads.first_read_file.input_fp);
			if(global_context->input_reads.is_paired_end_reads)
				fgets(read_line_buf2, 3000, global_context->input_reads.second_read_file.input_fp);

			flag = strtok_r(read_line_buf,"\t",&tok_tmp);
			if(!flag) break;

			flag = strtok_r(NULL,"\t",&tok_tmp);
			if(!flag) break;

			
			int flagi2 = 0;
			if(global_context->input_reads.is_paired_end_reads)
			{
				flag2 = strtok_r(read_line_buf2,"\t",&tok_tmp);
				flag2 = strtok_r(NULL,"\t",&tok_tmp);
				if(!flag2) break;
				flagi2 = atoi(flag2);
			}

			int flagi1 = atoi(flag);

			if((flagi1 & 0x100) == 0 && (flagi2&0x100) == 0)
			{
				read_position_1[processed_reads] = fhead_pos1;
				if(global_context->input_reads.is_paired_end_reads)
					read_position_2[processed_reads] = fhead_pos2;
				processed_reads++;
			}
			if(global_context->input_reads.is_paired_end_reads)
			{
				fgets(read_line_buf, 3000, global_context->input_reads.first_read_file.input_fp);
				fgets(read_line_buf2, 3000, global_context->input_reads.second_read_file.input_fp);
			}
		}
		//printf("PPPP=%llu\n", processed_reads);
		
	}
	else{
		while(1)
		{
			if(processed_reads >= chunk_reads || feof(global_context->input_reads.first_read_file.input_fp))
				break;

			read_position_1[processed_reads] = ftello(global_context->input_reads.first_read_file.input_fp);
			if(global_context->input_reads.is_paired_end_reads)
				read_position_2[processed_reads] = ftello(global_context->input_reads.second_read_file.input_fp);

			processed_reads++;

			geinput_jump_read(&global_context->input_reads.first_read_file);
			if(global_context->input_reads.is_paired_end_reads)
				geinput_jump_read(&global_context->input_reads.second_read_file);
		}
	}

	free(read_line_buf);
	free(read_line_buf2);

	int thread_no;
	for(thread_no = 0; thread_no < global_context->config.all_threads; thread_no++)
	{
		unsigned int my_start_read_no = processed_reads / global_context->config.all_threads * thread_no;
		unsigned int my_reads = (thread_no == global_context->config.all_threads-1)?(processed_reads - my_start_read_no):(processed_reads / global_context->config.all_threads);
		unsigned long long my_first_file_start = read_position_1[my_start_read_no]; 
		unsigned long long my_second_file_start = 0;
		if(global_context->input_reads.is_paired_end_reads)
			my_second_file_start = read_position_2[my_start_read_no];

		global_context -> input_reads.first_file_blocks[thread_no] = my_first_file_start;
		global_context -> input_reads.reads_in_blocks[thread_no] = my_reads;
		global_context -> input_reads.start_read_number_blocks[thread_no] = my_start_read_no;

		if(global_context->input_reads.is_paired_end_reads)
			global_context -> input_reads.second_file_blocks[thread_no] = my_second_file_start;
	}

	free(read_position_1);
	if(read_position_2)
		free(read_position_2);
	return processed_reads;
}

void locate_read_files(global_context_t * global_context, int type)
{
	if(type==SEEK_SET)
	{
		global_context -> current_circle_start_position_file1 = ftello(global_context -> input_reads.first_read_file.input_fp);
		if(global_context ->input_reads.is_paired_end_reads)
			global_context -> current_circle_start_position_file2 = ftello(global_context -> input_reads.second_read_file.input_fp);
	}
	else
	{
		global_context -> current_circle_end_position_file1 = ftello(global_context -> input_reads.first_read_file.input_fp);
		if(global_context ->input_reads.is_paired_end_reads)
			global_context -> current_circle_end_position_file2 = ftello(global_context -> input_reads.second_read_file.input_fp);
	
	}
}
void reward_read_files(global_context_t * global_context, int type)
{
	if(type==SEEK_SET)
	{
		fseeko(global_context -> input_reads.first_read_file.input_fp, global_context -> current_circle_start_position_file1, SEEK_SET);
		if(global_context ->input_reads.is_paired_end_reads)
			fseeko(global_context -> input_reads.second_read_file.input_fp, global_context -> current_circle_start_position_file2, SEEK_SET);
	}
	else
	{
		fseeko(global_context -> input_reads.first_read_file.input_fp, global_context -> current_circle_end_position_file1, SEEK_SET);
		if(global_context ->input_reads.is_paired_end_reads)
			fseeko(global_context -> input_reads.second_read_file.input_fp, global_context -> current_circle_end_position_file2, SEEK_SET);
	
	}
}


int read_chunk_circles(global_context_t *global_context)
{
	int block_no;

//	printf("GINP1 AT %llu\n", ftello(global_context -> input_reads.first_read_file.input_fp));
	
	while(1)
	{
		int ret;

		locate_read_files(global_context, SEEK_SET);
		if(global_context -> config.all_threads>1)
		{
			split_read_files(global_context);
			locate_read_files(global_context, SEEK_END);
			reward_read_files(global_context, SEEK_SET);
		}

		global_context -> current_index = (gehash_t*) malloc(sizeof(gehash_t));
		global_context -> current_value_index = (gene_value_index_t*) malloc(sizeof(gene_value_index_t));
		for(global_context->current_index_block_number = 0; global_context->current_index_block_number < global_context->index_block_number; global_context->current_index_block_number++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.tab", global_context->config.index_prefix, global_context->current_index_block_number,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			print_in_box(80,0,0, "Load the %d-th index block...",1+ global_context->current_index_block_number);


			if(gehash_load(global_context -> current_index, tmp_fname)) return -1;

			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, global_context->current_index_block_number, global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(global_context -> current_value_index, tmp_fname)) return -1;


			if(global_context->current_index_block_number ==0 && global_context -> all_processed_reads==0)
				global_context->align_start_time = miltime();

			ret = run_maybe_threads(global_context, STEP_VOTING);

			if(global_context -> config.all_threads<2 && global_context->current_index_block_number ==0)
				locate_read_files(global_context, SEEK_END);
			if(global_context->current_index_block_number < global_context->index_block_number -1)
				reward_read_files(global_context, SEEK_SET);

			gehash_destory_fast(global_context -> current_index);
			gvindex_destory(global_context -> current_value_index);
			if(ret) break;
			if(!global_context -> processed_reads_in_chunk) break;
		}

		free(global_context -> current_index);
		free(global_context -> current_value_index);


		//sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "%d reads have been processed in this chunk.", global_context -> processed_reads_in_chunk);


		// after the voting step, all subread index blocks are released and all base index blocks are loaded at once.

		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
		}

		if(!global_context -> processed_reads_in_chunk)
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			break;
	
		reward_read_files(global_context, SEEK_SET);
		ret = run_maybe_threads(global_context, STEP_ITERATION_ONE);
		ret = anti_supporting_read_scan(global_context);

		//HashTable * event_table = ((indel_context_t *)global_context -> module_contexts[MODULE_INDEL_ID])->event_entry_table;
		//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "There are %ld elements in the indel table before filtering.", event_table ->numOfElements);

		remove_neighbour(global_context);
		//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO, "There are only %ld elements in the indel table after filtering.", event_table ->numOfElements);

		reward_read_files(global_context, SEEK_SET);
		ret = ret || run_maybe_threads(global_context, STEP_ITERATION_TWO);

	//	printf("IBytes=%d+%d = %d\n", global_context -> all_value_indexes[0].start_base_offset, global_context -> all_value_indexes[0].values_bytes, global_context -> all_value_indexes[0].values [global_context -> all_value_indexes[0].values_bytes-1]);

		//gene_value_index_t * value_index = &global_context->all_value_indexes[0] ;
		
		//printf("=== I=%016llX   B=%016llX\n", (long long)value_index , (long long)value_index -> values);
		if(global_context -> config.is_third_iteration_running)
		{
			reward_read_files(global_context, SEEK_SET);
			ret = ret || do_iteration_three(global_context, NULL);
		}

		if(global_context -> config.report_sam_file)
		{
			reward_read_files(global_context, SEEK_SET);
			print_in_box(80, 0, 0, "%u %s were processed. Save the mapping results for them...", global_context ->processed_reads_in_chunk, global_context -> input_reads.is_paired_end_reads?"fragments":"reads");
			ret = ret || write_chunk_results(global_context);
			if('\r' == CORE_SOFT_BR_CHAR)
				sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO,"");
			
		}

		reward_read_files(global_context, SEEK_END);

		global_context -> all_processed_reads+= global_context ->processed_reads_in_chunk;

		if(ret) return ret;

		if(global_context -> processed_reads_in_chunk < global_context->config.reads_per_chunk)
			break; 
		else
			// base value indexes loaded in the last circle are not destroyed and are used in writting the indel VCF.
			// the indexes will be destroyed in destroy_global_context
			for(block_no = 0; block_no< global_context->index_block_number; block_no++)
				gvindex_destory(&global_context -> all_value_indexes[block_no]);

		
		clean_context_after_chunk(global_context);
	}

	// load all array index blocks at once.
	if(global_context -> config.is_third_iteration_running)
	{
		/*
		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
		{
			char tmp_fname[MAX_FILE_NAME_LENGTH];
			sprintf(tmp_fname, "%s.%02d.%c.array", global_context->config.index_prefix, block_no,  global_context->config.space_type == GENE_SPACE_COLOR?'c':'b');
			if(gvindex_load(&global_context -> all_value_indexes[block_no], tmp_fname)) return -1;
		}
		*/

		finalise_long_insertions(global_context);

		/*
		for(block_no = 0; block_no< global_context->index_block_number; block_no++)
			gvindex_destory(&global_context -> all_value_indexes[block_no]);
		*/
	}
	return 0;
}

void char_strftime(char * tbuf){
	time_t rawtime;
	struct tm * timeinfo;

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (tbuf,80,"%d-%b-%Y %H:%M:%S",timeinfo);

}

void print_subread_logo()
{
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m    _____ _    _ ____  _____  ______          _____  ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m =====      %c[0m%c[36m   / ____| |  | |  _ \\|  __ \\|  ____|   /\\   |  __ \\ ", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m   =====    %c[0m%c[36m  | (___ | |  | | |_) | |__) | |__     /  \\  | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m     ====   %c[0m%c[36m   \\___ \\| |  | |  _ <|  _  /|  __|   / /\\ \\ | |  | |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m       ==== %c[0m%c[36m   ____) | |__| | |_) | | \\ \\| |____ / ____ \\| |__| |", CHAR_ESC, CHAR_ESC, CHAR_ESC);
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %c[44;37m ========== %c[0m%c[36m  |_____/ \\____/|____/|_|  \\_\\______/_/    \\_\\_____/%c[0m", CHAR_ESC, CHAR_ESC, CHAR_ESC, CHAR_ESC);
	#ifdef MAKE_STANDALONE
	char * spaces = "";
	if(strlen(SUBREAD_VERSION) == 8) spaces = "";
	else if(strlen(SUBREAD_VERSION) == 5) spaces = "  ";
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"	%sv%s",spaces,SUBREAD_VERSION);
	#else
	sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_INFO ,"       %s",SUBREAD_VERSION);
	#endif
}



int print_configuration(global_context_t * context)
{
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
        print_subread_logo();
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");
        print_in_box(80, 1, 1, context->config.entry_program_name == CORE_PROGRAM_SUBJUNC?"subjunc setting":"subread-align setting");
        print_in_box(80, 0, 1, "");

        if(context->config.is_rna_seq_reads)
        {
                if(context->config.do_fusion_detection)
                {
                        print_in_box(80, 0, 0,         "          Function : Read alignment + Junction/Fusion detection%s", context->config.prefer_donor_receptor_junctions?" (RNA-Seq)":" (DNA-Seq)");
                }
                else
                        print_in_box(80, 0, 0,         "          Function : Read alignment + Junction detection (RNA-Seq)");
        }
        else
                print_in_box(80, 0, 0,         "          Function : Read alignment");
        print_in_box(80, 0, 0,         "           Threads : %d", context->config.all_threads);
        if( context->config.second_read_file[0])
        {
                print_in_box(80, 0, 0, "      Input file 1 : %s", context->config.first_read_file);
                print_in_box(80, 0, 0, "      Input file 2 : %s", context->config.second_read_file);
        }
        else
                print_in_box(80, 0, 0, "        Input file : %s%s", context->config.first_read_file, context->config.is_SAM_file_input?(context->config.is_BAM_input?" (BAM)":" (SAM)"):"");

        if(context->config.output_prefix [0])
                print_in_box(80, 0, 0, "       Output file : %s (%s)", context->config.output_prefix, context->config.is_BAM_output?"BAM":"SAM");
        else
                print_in_box(80, 0, 0, "     Output method : STDOUT (%s)" , context->config.is_BAM_output?"BAM":"SAM");

        print_in_box(80, 0, 0,         "        Index name : %s", context->config.index_prefix);
        print_in_box(80, 0, 0,         "      Phred offset : %d", (context->config.phred_score_format == FASTQ_PHRED33)?33:64);
        //print_in_box(80, 0, 0,         "        Space type : %s", (context->config.space_type == GENE_SPACE_COLOR)?"color-space":"base-space");
        print_in_box(80, 0, 1, "");
        if( context->config.second_read_file[0])
        {
                print_in_box(80, 0, 0, "   Min read1 votes : %d", context->config.minimum_subread_for_first_read);
                print_in_box(80, 0, 0, "   Min read2 votes : %d", context->config.minimum_subread_for_second_read);
                print_in_box(80, 0, 0, " Max fragment size : %d", context->config.maximum_pair_distance);
                print_in_box(80, 0, 0, " Min fragment size : %d", context->config.minimum_pair_distance);
                print_in_box(80, 0, 1, "");
        }
        else
                print_in_box(80, 0, 0, "         Min votes : %d", context->config.minimum_subread_for_first_read);

        print_in_box(80, 0, 0,         "        Max indels : %d", context->config.max_indel_length);
        print_in_box(80, 0, 0,         " # of Best mapping : %d", context->config.multi_best_reads);
        print_in_box(80, 0, 0,         "    Unique mapping : %s", context->config.report_multi_mapping_reads?"no":"yes");
        print_in_box(80, 0, 0,         "  Hamming distance : %s", context->config.use_hamming_distance_break_ties?"yes":"no");
        print_in_box(80, 0, 0,         "    Quality scores : %s", context->config.use_quality_score_break_ties?"yes":"no");

        if(context->config.max_insertion_at_junctions)
                print_in_box(80, 0, 0,         "Insertions at junc : %d", context->config.max_insertion_at_junctions);

        if(context->config.read_group_id[0])
                print_in_box(80, 0, 0, "   Read group name : %s", context->config.read_group_id);

        print_in_box(80, 0, 1, "");
        print_in_box(80, 2, 1, "http://subread.sourceforge.net/");
        sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"");

        if(!context->config.first_read_file[0])
        {
                sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify at least one input file in the FASTQ/FASTA/PLAIN format using the '-r' option.\n");
                return -1;
        }

        if(0 && !context->config.output_prefix[0])
        {
                sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the path of output using the '-o' option.\n");
                return -1;
        }

        if(!context->config.index_prefix[0])
        {
                sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"You have to specify the prefix of the index files using the '-i' option.\n");
                return -1;
        }
        char tbuf[90];
        char_strftime(tbuf);

        print_in_box(80,1,1,"Running (%s)", tbuf);
        print_in_box(80,0,1,"");


        return 0;
}



int init_paired_votes(global_context_t *context)
{

	if(context -> config.is_rna_seq_reads)
		context -> chunk_subjunc_records = malloc(sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);
	else	context -> chunk_subjunc_records = NULL;
	context -> chunk_alignment_records = malloc(sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);


	if(!context -> chunk_alignment_records)
	{	
		return 1;
	}

	context -> big_margin_record = malloc( sizeof(*context -> big_margin_record) * (context->input_reads.is_paired_end_reads?2:1) * context -> config.big_margin_record_size * context ->config.reads_per_chunk);

	memset(context ->big_margin_record  , 0 , sizeof(*context -> big_margin_record) *context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context -> config.big_margin_record_size);
	memset(context ->chunk_alignment_records , 0 , sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);

		//fprintf(stderr, "MALLOC=%llu = %d * %d * %d \n", sizeof(alignment_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads, sizeof(alignment_result_t), context ->config.reads_per_chunk, context->config.multi_best_reads);
		//sleep(10000);

	if(context -> chunk_subjunc_records)
		memset(context ->chunk_subjunc_records , 0 , sizeof(subjunc_result_t) * context ->config.reads_per_chunk * (context->input_reads.is_paired_end_reads?2:1) * context->config.multi_best_reads);


	return 0;
}

void write_sam_headers(global_context_t * context)
{
	if(context -> config.is_BAM_output)
	{
		SamBam_writer_add_header(context -> output_bam_writer,"@HD\tVN:1.0\tSO:unsorted", 0);
		int xk1;
		unsigned int last_offset = 0;
		char obuf[300];
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			SamBam_writer_add_chromosome(context -> output_bam_writer, context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, context->chromosome_table.read_offsets[xk1] - last_offset+16, 1);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}


		if(context->config.read_group_id[0])
		{
			snprintf(obuf,299, "@RG\tID:%s%s",context->config.read_group_id, context->config.read_group_txt);
			SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
		}
		snprintf(obuf,299, "@PG\tID:subread\tPN:subread\tVN:%s", SUBREAD_VERSION);
		SamBam_writer_add_header(context -> output_bam_writer,obuf, 0);
	}
	else
	{
		sambamout_fprintf(context -> output_sam_fp, "@HD\tVN:1.0\tSO:unsorted\n");
		int xk1;
		unsigned int last_offset = 0;
		for(xk1=0; xk1< context->chromosome_table.total_offsets; xk1++)
		{
			sambamout_fprintf(context -> output_sam_fp, "@SQ\tSN:%s\tLN:%u\n", context->chromosome_table.read_names+ xk1 * MAX_CHROMOSOME_NAME_LEN, context->chromosome_table.read_offsets[xk1] - last_offset+16);
			last_offset = context->chromosome_table.read_offsets[xk1];
		}

		if(context->config.read_group_id[0])
			sambamout_fprintf(context -> output_sam_fp, "@RG\tID:%s%s\n",context->config.read_group_id, context->config.read_group_txt);
		sambamout_fprintf(context -> output_sam_fp, "@PG\tID:subread\tPN:subread\tVN:%s\n", SUBREAD_VERSION);
		
	}
}

int load_global_context(global_context_t * context)
{
	char tmp_fname [MAX_FILE_NAME_LENGTH];
	int min_phred_score = -1 , max_phred_score = -1;
	context -> is_phred_warning = 0; 
	

	if(core_geinput_open(context, &context->input_reads.first_read_file, 1,1))
	{
		//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' as input. Please check if it exists, you have the permission to read it, and it is in the correct format.\n", context->config.first_read_file);
		return -1;
	}

	context->config.space_type = context->input_reads.first_read_file.space_type;
	print_in_box(89,0,0,"The input file contains %c[36m%s%c[0m space reads.", CHAR_ESC, context->config.space_type == GENE_SPACE_COLOR?"color":"base", CHAR_ESC);
	if(context->config.space_type == GENE_SPACE_COLOR && context->config.is_BAM_output && !context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The color-space bases will be converted into base space in the BAM output.");
		context->config.convert_color_to_base=1;
	}
	else if(context->config.space_type == GENE_SPACE_BASE && context->config.convert_color_to_base)
	{
		print_in_box(80,0,0,"The reads will not be converted into base space.");
		context->config.convert_color_to_base=0;
	}

	if(context->input_reads.is_paired_end_reads)
	{
		if(core_geinput_open(context, &context->input_reads.second_read_file, 2,1))
		{
			//sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' as input. Please check if it exists, you have the permission to read it, and it is in the correct format.\n", context->config.second_read_file);
			return -1;
		}
	}


	if(context->input_reads.is_paired_end_reads)
		context->config.reads_per_chunk = 7*1024*1024*min(40,max(0.01,context->config.memory_use_multiplex));
	else
		context->config.reads_per_chunk = 14*1024*1024*min(40,max(0.01,context->config.memory_use_multiplex));

	struct stat ginp1_stat;
	stat(context->config.first_read_file , &ginp1_stat);
	context->input_reads.first_read_file_size = ginp1_stat.st_size;


	context -> input_reads.avg_read_length = guess_reads_density_format(context->config.first_read_file , context->config.is_SAM_file_input?1:0, &min_phred_score, &max_phred_score);
	if(context -> input_reads.avg_read_length<0 )context -> input_reads.avg_read_length = 250;
//	SUBREADprintf("QR=[%d,%d]; ALEN=%f\n",  min_phred_score, max_phred_score, context -> input_reads.avg_read_length);
	if(max_phred_score>=0)
	{
		if((context->config.phred_score_format == FASTQ_PHRED64 && min_phred_score < 65) || (context->config.phred_score_format == FASTQ_PHRED33 && max_phred_score > 33+50))
		{
			print_in_box(80,0,0, "WARNING The specified phred-score offset (%d) seems to be incorrect.", context->config.phred_score_format == FASTQ_PHRED33?33:64);
			print_in_box(80,0,0, "        The observed phred-score range is [%d,%d].", min_phred_score, max_phred_score);
			print_in_box(80,0,0, "");
			context -> is_phred_warning = 1;
		}
	}


	if(context->config.report_sam_file && context -> config.output_prefix[0])
	{
		// ====== open output files ======
		// Only the sam file is opened here; other files like bed, indel and etc are opened in init_modules()
		sprintf(tmp_fname,"%s", context->config.output_prefix);

		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer , tmp_fname);
			context -> output_sam_fp = NULL;
		}
		else
		{
			context -> output_sam_fp = f_subr_open(tmp_fname,"wb");
			context -> output_bam_writer = NULL;
		}
		if((!context -> output_bam_writer) && (!context->output_sam_fp))
		{
			sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable to open '%s' to write. Please check if the path exists and you have the permission to create/write this file", tmp_fname);
			return -1;
		}
	}
	else
	{
		if(context -> config.is_BAM_output)
		{
			context -> output_bam_writer = malloc(sizeof(SamBam_Writer));
			SamBam_writer_create(context -> output_bam_writer ,NULL);
		}
		context->output_sam_fp = NULL;
	}
	
	// ====== check index files, count blocks and load chro table ======
	sprintf(tmp_fname, "%s.reads", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Unable top open index '%s'. Please make sure that the correct prefix is specified and you have the permission to read these files. For example, if there are files '/opt/my_index.reads', '/opt/my_index.files' and etc, the index prefix should be specified as '/opt/my_index' without any suffix. \n", context->config.index_prefix);
		return -1;
	}

	if(context->config.space_type == GENE_SPACE_COLOR)
		sprintf(tmp_fname, "%s.00.c.tab", context->config.index_prefix);
	else
		sprintf(tmp_fname, "%s.00.b.tab", context->config.index_prefix);
	if(!does_file_exist(tmp_fname))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Your reads are in the %s space but the index was not built in the same space. Unable to precess the reads.\n", context->config.space_type == GENE_SPACE_COLOR?"color":"base");
		return -1;
	}


	context->index_block_number = 0; 
	while(1)
	{
		sprintf(tmp_fname, "%s.%02d.%c.tab", context->config.index_prefix, context->index_block_number, context->config.space_type == GENE_SPACE_COLOR?'c':'b');
		if(!does_file_exist(tmp_fname))break;
		context->index_block_number ++;
		if(context->index_block_number>=2 && context->config.max_indel_length > 16)
		{
			print_in_box(80,0,0,"ERROR You cannot use multi-block index for very-long indel detection!");
			print_in_box(80,0,0,"Please set the maximum indel length <= 16.");
			return -1;
		}
	}

	context->current_index_block_number = 0;
	load_offsets(&context->chromosome_table, context->config.index_prefix);

	if(context->config.report_sam_file)
		write_sam_headers(context);

	// ====== init other variables ======
	if(context -> config.all_threads>1)
		subread_init_lock(&context -> thread_initial_lock);

	if(init_paired_votes(context))
	{
		sublog_printf(SUBLOG_STAGE_RELEASED, SUBLOG_LEVEL_ERROR,"Cannot initialise the voting space. You need at least 2GB of empty physical memory to run this program.\n");
		return 1;
	}
	context->all_processed_reads = 0;
	context->all_mapped_reads = 0;
	context->all_correct_PE_reads = 0;
	context->all_junctions = 0;
	context->all_fusions = 0;
	context->all_indels = 0;
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "load_global_context: finished");

	memset( context->all_value_indexes , 0 , 100 * sizeof(gene_value_index_t));

	return 0;
}

int init_modules(global_context_t * context)
{
	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: begin");
	int ret = init_indel_tables(context);
	if(context->config.is_rna_seq_reads || context->config.do_fusion_detection)
		ret = ret || init_junction_tables(context);

	sublog_printf(SUBLOG_STAGE_DEV1, SUBLOG_LEVEL_DEBUG, "init_modules: finished: %d",ret);
	return ret;
}

int destroy_modules(global_context_t * context)
{
	destroy_indel_module(context);
	if(context->config.is_rna_seq_reads)
		destroy_junction_tables(context);
	return 0;
}

int destroy_global_context(global_context_t * context)
{
	int xk1, block_no;

	for(block_no = 0; block_no< context->index_block_number; block_no++)
		gvindex_destory(&context -> all_value_indexes[block_no]);

	if(context->output_sam_fp)
		fclose(context->output_sam_fp);
	if(context->output_bam_writer)
	{
		SamBam_writer_close(context->output_bam_writer);
		free(context->output_bam_writer);
		context->output_bam_writer=NULL;
	}
	free(context->chunk_alignment_records);
	free(context->big_margin_record);
	if(context->chunk_subjunc_records)
		free(context->chunk_subjunc_records);
	
	for(xk1=0; xk1<5; xk1++)
		if(context->module_contexts[xk1])free(context->module_contexts[xk1]);
	geinput_close(&context -> input_reads.first_read_file);
	if(context->input_reads.is_paired_end_reads) geinput_close(&context -> input_reads.second_read_file);
	destroy_offsets(&context->chromosome_table);

	
	if((context -> will_remove_input_file & 1) && (memcmp(context ->config.first_read_file, "./core-temp", 11) == 0)) unlink(context ->config.first_read_file);
	if((context -> will_remove_input_file & 2) && (memcmp(context ->config.second_read_file, "./core-temp", 11) == 0)) unlink(context ->config.second_read_file);
	
	return 0;
}


int write_bincigar_part(char * bincigar, int chropt, unsigned int optlen, int bincigar_len)
{
	int binopt, binbytes, x1;

	if(optlen<1) return -1;

	if(optlen < 4) binbytes=1;
	else if(optlen < 1024) binbytes=2;
	else if(optlen < 262144) binbytes=3;
	else if(optlen < 67108864) binbytes=4;
	else binbytes=5;

	if(bincigar_len<binbytes) return -1; 

	switch(chropt)
	{
		case 'S':
			binopt = CORE_CIGAR_OPT_S;
			break;
		case 'M':
			binopt = CORE_CIGAR_OPT_M;
			break;
		case 'I':
			binopt = CORE_CIGAR_OPT_I;
			break;
		case 'D':
			binopt = CORE_CIGAR_OPT_D;
			break;
		case 'N':
			binopt = CORE_CIGAR_OPT_N;
			break;
		case 'n':
			binopt = CORE_CIGAR_OPT_NN;
			break;
		case 'B':
			binopt = CORE_CIGAR_OPT_B;
			break;
		case 'b':
			binopt = CORE_CIGAR_OPT_BB;
			break;
		default:
			return -1;
	}

	bincigar[0]=binopt | (binbytes << 3) | ((optlen & 3)<< 6);
	optlen >>= 2;
	for(x1=1;x1<binbytes; x1++)
	{
		bincigar[x1] = optlen&0xff;
		optlen>>=8;
	}

	return binbytes;
}

// function returns the actual length of bincigar, or -1 if anything is wrong, e.g., bincigar_len is too short or unrecognized operations.
int cigar2bincigar(char *cigar, char *bincigar, int bincigar_len)
{
	int xk1=0;
	unsigned int tmpv=0, bincigar_cursor=0;
	while(1)
	{
		int nch = cigar[xk1];
		if(!nch) break;
		xk1++;

		if(isdigit(nch)) tmpv=tmpv*10+(nch-'0');
		else
		{
			int bincigar_sec_len = write_bincigar_part(bincigar+bincigar_cursor, nch, tmpv, bincigar_len-bincigar_cursor);
			if(bincigar_sec_len<0){
				bincigar[0]=0;
				return -1;
			}
			bincigar_cursor += bincigar_sec_len;
			tmpv=0;
		}
	}

	if(bincigar_cursor<bincigar_len) bincigar[bincigar_cursor] = 0;

	//printf("%s : BL=%d\n", cigar, bincigar_cursor);

	return bincigar_cursor;
}


int write_cigar_part(char *bincigar, char *cigar, int cigar_len , int * bincigar_move)
{
	int binbytes, x1, binopt, charopt;
	unsigned int tmpv = 0;
	char sec_buf[13];

	binbytes = 7& (bincigar[0] >> 3);
	binopt = 7 & bincigar[0];

	switch(binopt)
	{
		case CORE_CIGAR_OPT_D:
			charopt='D'; 
			break;
		case CORE_CIGAR_OPT_I:
			charopt='I'; 
			break;
		case CORE_CIGAR_OPT_M:
			charopt='M'; 
			break;
		case CORE_CIGAR_OPT_S:
			charopt='S'; 
			break;
		case CORE_CIGAR_OPT_B:
			charopt='B'; 
			break;
		case CORE_CIGAR_OPT_BB:
			charopt='b'; 
			break;
		case CORE_CIGAR_OPT_N:
			charopt='N'; 
			break;
		default:
			charopt='n'; 
			break;
	}

	tmpv = (bincigar[0]>>6) & 3;
	for(x1 = 1; x1 < binbytes; x1++)
	{
		unsigned int dtmpv = 0xff & bincigar[x1];
		dtmpv <<= (x1*8 - 6);
		tmpv += dtmpv; 
	}

	int added_len = sprintf(sec_buf, "%u%c", tmpv, charopt);
	if(added_len > cigar_len)
		return -1;
	memcpy(cigar, sec_buf, added_len);
	(*bincigar_move) = binbytes;

	return added_len;
}

int bincigar2cigar(char * cigar, int cigar_len, char * bincigar, int bincigar_max_len, int read_len)
{
	int cigar_cursor = 0, bincigar_cursor = 0;
	while(1)
	{
		int bincigar_move = 0;
		int cigar_sec_len = write_cigar_part(bincigar + bincigar_cursor, cigar+cigar_cursor, cigar_len-cigar_cursor-1, &bincigar_move);
		if(cigar_sec_len<0){
			sprintf(cigar,"%dM", read_len);
			return -1;
		}
		//printf("NPC=%s\n", cigar);
		cigar_cursor += cigar_sec_len;
		bincigar_cursor += bincigar_move;
		if(bincigar_cursor>=bincigar_max_len) break;
		if(bincigar[bincigar_cursor] == 0) break;
	}
	cigar[cigar_cursor] = 0;

	return cigar_cursor;
}

int term_strncpy(char * dst, char * src, int max_dst_mem)
{
	int i;

	for(i=0; i<max_dst_mem; i++)
	{
		if(!src[i]) break;
		dst[i]=src[i];
		if(i == max_dst_mem-1)
			SUBREADprintf("String out of memory limit: '%s'\n", src);
	}
	if(i >= max_dst_mem) i = max_dst_mem-1;
	dst[i] = 0;

	return 0;
}


// This assumes the first part of Cigar has differet strandness to the main part of the cigar.
// Pos is the LAST WANTED BASE location before the first strand jump (split by 'b' or 'n').
// The first base in the read actually has a larger coordinate than Pos. 
// new_cigar has to be at least 100 bytes.
unsigned int reverse_cigar(unsigned int pos, char * cigar, char * new_cigar)
{
	int cigar_cursor = 0;
	new_cigar[0]=0;
	unsigned int tmpi=0;
	int last_piece_end = 0;
	int last_sec_start = 0;
	unsigned int chro_pos = pos, this_section_start = pos, ret = pos;
	int is_positive_dir = 0;
	int read_cursor = 0;
	int section_no = 0;

	for(cigar_cursor = 0 ;  ; cigar_cursor++)
	{
		if( cigar [cigar_cursor] == 'n' ||  cigar [cigar_cursor] == 'b' ||  cigar [cigar_cursor] == 0)
		{
			int xk1, jmlen=0, nclen=strlen(new_cigar);
			char jump_mode [13];

			if(cigar [cigar_cursor] !=0)
			{
				sprintf(jump_mode, "%u%c", tmpi,  cigar [cigar_cursor] == 'b'?'n':'b');
				jmlen = strlen(jump_mode);
			}

			for(xk1=nclen-1;xk1>=0; xk1--)
				new_cigar[ xk1 +  last_piece_end + jmlen - last_sec_start ] = new_cigar[ xk1 ];
			new_cigar [nclen + jmlen + last_piece_end - last_sec_start ] = 0;

			memcpy(new_cigar , jump_mode, jmlen);
			memcpy(new_cigar + jmlen , cigar + last_sec_start, last_piece_end - last_sec_start);

			last_sec_start = cigar_cursor+1;

			if(is_positive_dir && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos -= tmpi - read_cursor - 1;
				else	chro_pos += tmpi - read_cursor - 1;
			}
			if((!is_positive_dir) && cigar [cigar_cursor] !=0)
			{
				if(cigar [cigar_cursor] == 'b') chro_pos = this_section_start - tmpi - read_cursor - 1;
				else	chro_pos = this_section_start + tmpi - read_cursor - 1;
			}

			this_section_start = chro_pos;

			if(section_no == 0)
				ret = chro_pos;

			is_positive_dir = ! is_positive_dir;
			section_no++;
			tmpi=0;
		}
		else if(isalpha(cigar [cigar_cursor]))
		{
			if(cigar [cigar_cursor]=='M' || cigar [cigar_cursor] == 'S')
				read_cursor += tmpi;
			tmpi=0;
			last_piece_end = cigar_cursor+1;
		}
		else tmpi = tmpi*10 + (cigar [cigar_cursor] - '0');

		if(cigar [cigar_cursor] == 0)break;
	}

	//printf("REV CIGAR: %s  =>  %s\n", cigar, new_cigar);
	return ret;
}

int chimeric_cigar_parts(global_context_t * global_context, unsigned int sel_pos, int is_first_section_negative_strand, int is_first_section_reversed, char * in_cigar, unsigned int * out_poses, char ** out_cigars, char * out_strands, int read_len, short * perfect_lens)
{
	unsigned int current_perfect_map_start = sel_pos;
	int current_perfect_section_no = 0;
	int current_perfect_cursor = sel_pos;
	int is_reversed = is_first_section_reversed;
	int is_negative = is_first_section_negative_strand;
	int read_cursor = 0;
	int out_cigar_writer_ptr = 0;
	unsigned int tmpi = 0;

	short perfect_len = 0;

	int cigar_cursor;

	out_poses[0] = current_perfect_map_start;
	out_strands[0] = is_negative?'-':'+';
	char main_piece_strand = (is_first_section_negative_strand == is_first_section_reversed)?'+':'-';

	for(cigar_cursor=0;;cigar_cursor++)
	{
		char ncch = in_cigar[cigar_cursor];
		int is_chimeric_section_end = 0;

		if(!ncch){
			perfect_lens [current_perfect_section_no] = perfect_len ;
			current_perfect_section_no++;
			break;
		}

		if(toupper(ncch)=='N'||toupper(ncch)=='B')
		{

			unsigned int jummped_location;
			int is_chro_jump = 0, is_long_jump = 0;

			if(is_reversed)
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_map_start - 1 + tmpi;
				else
					jummped_location = current_perfect_map_start - 1 - tmpi;
			}
			else
			{
				if(toupper(ncch)=='N')
					jummped_location = current_perfect_cursor + tmpi;
				else
					jummped_location = current_perfect_cursor - tmpi;

			}

			if(ncch == 'N')
			{
				char * curr_chr, * new_chr;
				unsigned int curr_offset, new_offset;
				locate_gene_position_max(current_perfect_cursor, &global_context -> chromosome_table, & curr_chr, & curr_offset, 1);
				locate_gene_position_max(jummped_location      , &global_context -> chromosome_table, &  new_chr, &  new_offset, 1);
				assert(curr_chr);
				assert(new_chr);
				is_chro_jump = (curr_chr != new_chr);

				long long int dist = current_perfect_cursor;
				dist -= jummped_location;
				if(abs(dist) >= 134217728)
					is_long_jump = 1;
				// A long jump is the jump longer than 2^27.
				// Picard does not like it!!
			}

			if(is_chro_jump || islower(ncch) || ncch == 'B' || is_long_jump)
			{
				current_perfect_cursor = jummped_location;

				if(islower(ncch)){
					is_reversed = !is_reversed;
					is_negative = !is_negative;
				}

				current_perfect_map_start = current_perfect_cursor;
				tmpi = 0;
				if(read_cursor<read_len)
					sprintf(out_cigars[current_perfect_section_no] + out_cigar_writer_ptr,"%dS", read_len - read_cursor);

				perfect_lens [current_perfect_section_no] = perfect_len ;
				perfect_len = 0;

				current_perfect_section_no++;
				if(current_perfect_section_no>CIGAR_PERFECT_SECTIONS)break;

				out_poses[current_perfect_section_no] = current_perfect_map_start - read_cursor;
				out_strands[current_perfect_section_no] = is_negative?'-':'+';
				out_cigar_writer_ptr = sprintf(out_cigars[current_perfect_section_no],"%dS", read_cursor);
				is_chimeric_section_end  = 1;
			}
		}

		if(!is_chimeric_section_end)
		{
			if(isalpha(ncch))
			{
				out_cigar_writer_ptr+=sprintf(out_cigars[current_perfect_section_no]+out_cigar_writer_ptr, "%u%c", tmpi, ncch);
			}
			if(ncch == 'M'|| ncch == 'S')
			{
				read_cursor += tmpi;
				if(ncch == 'M')
					perfect_len += tmpi;
				if(!is_reversed)
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'D' || ncch == 'N')
			{
				if(!is_reversed)
					current_perfect_cursor += tmpi;
				tmpi = 0;
			}
			else if(ncch == 'I')
			{
				read_cursor += tmpi;
				tmpi = 0;
			}
			else if(isdigit(ncch))
				tmpi = tmpi*10+(ncch-'0');
		}
	}

	int xk1 = 0, best_match = -9999;

	for(xk1=0; xk1<current_perfect_section_no;xk1++)
	{
		if(best_match < 0 || (main_piece_strand == out_strands[xk1] && perfect_lens[xk1]>perfect_lens[best_match]))
			best_match = xk1;
	}

	if(best_match>0)
	{
		unsigned int tmpv;
		char cigar_sec[100];
		tmpv = out_poses[0];
		out_poses[0]=out_poses[best_match];
		out_poses[best_match] = tmpv;

		tmpv = out_strands[0];
		out_strands[0] = out_strands[best_match];
		out_strands[best_match] = tmpv;

		strcpy(cigar_sec, out_cigars[0]);
		strcpy(out_cigars[0], out_cigars[best_match]);
		strcpy(out_cigars[best_match] , cigar_sec);

		tmpv = perfect_lens[0];
		perfect_lens[0] = perfect_lens[best_match];
		perfect_lens[best_match] = tmpv;
	}

	return current_perfect_section_no;
}

void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r));

void quick_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	quick_sort_run(arr, 0, arr_size-1, compare, exchange);
}
 
 
void quick_sort_run(void * arr, int spot_low,int spot_high, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r))
{
	int pivot,j,i;

	if(spot_high-spot_low<1) return;
	pivot = spot_low;
	i = spot_low;
	j = spot_high;

	while(i<=j)
	{
		if(compare(arr, i, pivot) <0)
		{
			i++;
			continue;
		}

		if(compare(arr, j, pivot)>0)
		{
			j--;
			continue;
		}

		if(i!=j)
			exchange(arr,i,j);
		i++;
		j--;
	}

	quick_sort_run(arr, spot_low, j, compare, exchange);
	quick_sort_run(arr, i, spot_high, compare, exchange);
	
}



void merge_sort_run(void * arr, int start, int items, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	if(items > 11)
	{
		int half_point = items/2;
		merge_sort_run(arr, start, half_point, compare, exchange, merge);
		merge_sort_run(arr, start + half_point, items - half_point, compare, exchange, merge);
		merge(arr, start, half_point, items - half_point);
	}
	else
	{
		int i, j;
		for(i=start; i< start + items - 1; i++)
		{
			int min_j = i;
			for(j=i + 1; j< start + items; j++)
			{
				if(compare(arr, min_j, j) > 0)
					min_j = j;
			}
			if(i!=min_j)
				exchange(arr, i, min_j);
		}
	}
}
void merge_sort(void * arr, int arr_size, int compare (void * arr, int l, int r), void exchange(void * arr, int l, int r), void merge(void * arr, int start, int items, int items2))
{
	merge_sort_run(arr, 0, arr_size, compare, exchange, merge);
}
