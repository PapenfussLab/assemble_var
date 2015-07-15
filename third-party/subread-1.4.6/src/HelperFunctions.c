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
  
  
#include <ctype.h>
#include <string.h>

#include "subread.h"
#include "HelperFunctions.h"

int RSubread_parse_CIGAR_string(const char * CIGAR_Str, unsigned int * Staring_Points, unsigned short * Section_Length)
{
	unsigned int tmp_int=0;
	int cigar_cursor=0;
	unsigned short read_cursor=0;
	unsigned int chromosome_cursor=0;
	int ret=0;

	for(cigar_cursor=0; ; cigar_cursor++)
	{
		char ch = CIGAR_Str[cigar_cursor];

		if(ch >='0' && ch <= '9')
		{
			tmp_int=tmp_int*10+(ch - '0');
		}
		else
		{
			if(ch == 'M' || ch == 'D')
			{
				read_cursor += tmp_int;
				chromosome_cursor += tmp_int;
			}
			else if(ch == 'N' || ch == 0)
			{
				if(ret <6)
				{
					if(read_cursor>0)
					{
						Staring_Points[ret] = chromosome_cursor - read_cursor;
						Section_Length[ret] = read_cursor;
						ret ++;
					}
				}
				read_cursor = 0;

				if(ch == 'N') chromosome_cursor += tmp_int;
				else break;
			}
			//printf("C=%c, TV=%d, CC=%d, RC=%d\n", ch, tmp_int, chromosome_cursor, read_cursor);
			tmp_int = 0;
		}
		if(cigar_cursor>100) return -1;
	}

	return ret;
}

void display_sections(char * CIGAR_Str)
{
	unsigned int Staring_Points[6];
	unsigned short Section_Length[6];
	int retv = RSubread_parse_CIGAR_string(CIGAR_Str, Staring_Points, Section_Length);

	int x1;
	SUBREADprintf("Cigar=%s ; Sections=%d\n", CIGAR_Str, retv);
	for(x1=0; x1<retv; x1++)
	{
		SUBREADprintf("   Section #%d: offset=%u  length=%u\n",x1, Staring_Points[x1], Section_Length[x1]);
	}
	SUBREADprintf("\n");
	
}


#define GECV_STATE_BEFORE 10
#define GECV_STATE_NAME 20
#define GECV_STATE_GAP 30
#define GECV_STATE_VALUE 40
#define GECV_STATE_QVALUE 50
#define GECV_STATE_QV_END 60
#define GECV_STATE_ERROR 9999

int GTF_extra_column_istoken_chr(char c)
{
	return (isalpha(c)||isdigit(c)||c=='_');
}

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size)
{
	int state = GECV_STATE_BEFORE;
	int col_cursor = 0, is_correct_name=0;
	char name_buffer[200];
	int name_cursor = 0, value_cursor=-1;

	while(1)
	{
		if(name_cursor>190) return -1;
		char nch = Extra_Col[col_cursor];
		if(nch == '\n' || nch == '\r') nch = 0;
		if(state == GECV_STATE_BEFORE)
		{
			if(GTF_extra_column_istoken_chr(nch))
			{
				name_buffer[0] = nch;
				name_cursor = 1;
				state = GECV_STATE_NAME;
			}
			else if(nch != ' ' && nch != 0)
			{
				state = GECV_STATE_ERROR;
			}
		}
		else if(state == GECV_STATE_NAME)
		{
			if(nch == ' ' || nch == '=')
			{
				state = GECV_STATE_GAP;
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				//printf("CORR=%d : '%s'\n", is_correct_name, name_buffer);
			}
			else if(nch == '"')
			{
				name_buffer[name_cursor] = 0;
				is_correct_name = (strcmp(name_buffer , Target_Name) == 0);
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(GTF_extra_column_istoken_chr(nch))
				name_buffer[name_cursor++] = nch;
			else
			{
				state = GECV_STATE_ERROR;
				//printf("ERR2  : '%c'\n", nch);
			}
			
		}
		else if(state == GECV_STATE_GAP)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QVALUE;
				if(is_correct_name)
					value_cursor = 0;
			}
			else if(nch != '=' && isgraph(nch))
			{
				state = GECV_STATE_VALUE;
				if(is_correct_name)
				{
					Target_Value[0]=nch;
					value_cursor = 1;
				}
			}
			else if(nch != ' ' && nch != '=')
				state = GECV_STATE_ERROR;
		}
		else if(state == GECV_STATE_VALUE)
		{
			if(nch == ';' || nch == 0)
			{
				state = GECV_STATE_BEFORE;
				if(is_correct_name)
				{
					Target_Value[value_cursor] = 0;
				}
				is_correct_name = 0;
			}
			else{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
					Target_Value[value_cursor++] = nch;
			}
		}
		else if(state == GECV_STATE_QVALUE)
		{
			if(nch == '"')
			{
				state = GECV_STATE_QV_END;
				if(is_correct_name)
					Target_Value[value_cursor] = 0;
				is_correct_name = 0;
			}
			else
			{
				if(value_cursor < TargVal_Size-1 && is_correct_name)
				{
					if(nch !=' ' || value_cursor>0)
						Target_Value[value_cursor++] = nch;
				}
			}
		}
		else if(state == GECV_STATE_QV_END)
		{
			if(nch == ';' || nch == 0)
				state = GECV_STATE_BEFORE;
			else if(nch != ' ')
				state = GECV_STATE_ERROR;
				
		}

		if (GECV_STATE_ERROR == state){
			Target_Value[0]=0;
			return -1;
		}
		if (nch == 0)
		{
			if(state == GECV_STATE_BEFORE && value_cursor>0)
			{
				int x1;
				for(x1 = value_cursor-1; x1>=0; x1--)
				{
					if(Target_Value[x1] == ' '){
						value_cursor --;
						Target_Value[x1]=0;
					}
					else break;
				}

				if(value_cursor>0)
					return value_cursor;
			}
			Target_Value[0]=0;
			return -1;
		}
		col_cursor++;
	}
}


void hpl_test2_func()
{
	char * extra_column = " gene_id \"PC4-013  \"; 013=ABCD  ; PC4 =  CCXX  ";
	char * col_name = "gene_id";
	char col_val[100];

	int col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	extra_column = "gene_id =   \"PC4-013  ;=\"  ;013 = AXXD ; PC4=x";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);




	extra_column = " gene_id\"  PC4-013  ;=  \"; XXX='123' ;013 :ABCD  ; PC4 =  CCXX=  ;";
	col_name = "013";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "XXX";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);


	col_name = "PC4";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);

	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);



	extra_column = "gene_id \"653635\"; transcript_id \"TR:653635\";";
	col_name = "gene_id";
	col_len = GTF_extra_column_value(extra_column, col_name, col_val, 100);
	SUBREADprintf("LEN=%d; KEY='%s'; VAL=\"%s\"\n", col_len, col_name, col_val);





}

void hpl_test1_func()
{
	display_sections("");
	display_sections("*");
	display_sections("5S10M2D10M800N12M3I12M450N12M12D99M6S");
	display_sections("110M2I10M800N32M3I12M6S");
	display_sections("200S110M2I10M800N32M3I12M200N40M");
	display_sections("3M1663N61M1045N36M3D20M66N10M2D10M77N3M1663N61M1045N36M3D20M66N103M1663N61M1045N36M3D20M66N9M");
}

#ifdef RSUBREAD_TEST_HELPER_FUNCTIONS
void main()
#else
void testi_helper_1_main()
#endif
{
	hpl_test2_func();
}

char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep
    int len_with; // length of with
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    if (!orig)
        return NULL;
    if (!rep)
        rep = "";
    len_rep = strlen(rep);
    if (!with)
        with = "";
    len_with = strlen(with);

    ins = orig;
    for (count = 0; NULL != (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }
    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}



// rule: the string is ABC123XXXXXX...
// This is the priroity:
// First, compare the letters part.
// Second, compare the pure numeric part.
// Third, compare the remainder.
int strcmp_number(char * s1, char * s2)
{
	int x1 = 0;
	int ret = 0;

	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];

		if((!nch1) || !nch2){return nch2?1:(nch1?(-1):0);}
		if(isdigit(nch1) && isdigit(nch2))break;

		ret = nch1 - nch2;
		if(ret) return ret;
		x1++;
	}

	int v1 = 0, v2 = 0;
	while(1)
	{
		char nch1 = s1[x1];
		char nch2 = s2[x1];
		if((!nch1) || !nch2){
			if(nch1 || nch2)
				return nch2?(-1):1;
			break;
		}
		int is_chr1_digit = isdigit(nch1);
		int is_chr2_digit = isdigit(nch2);

		if(is_chr1_digit || is_chr2_digit)
		{
			if(is_chr1_digit && is_chr2_digit)
			{
				v1 = v1*10+(nch1-'0');
				v2 = v2*10+(nch2-'0');
			}
			else
			{
				ret = nch1 - nch2;
				return ret;
			}
		}
		else break;
		x1++;
	}

	if(v1==v2)
		return strcmp(s1+x1, s2+x1);	
	else
	{
		ret = v1 - v2;
		return ret;
	}
}
