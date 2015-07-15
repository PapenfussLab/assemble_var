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
  
  
#ifndef __HELPER_FUNCTIONS_H_
#define __HELPER_FUNCTIONS_H_


// This function parses CIGAR_Str and extract the relative starting points and lengths of all sections (i.e., the sections of read that are separated by 'N').
// CIGAR_Str is a CIGAR string containing 'S', 'M', 'I', 'D' and 'N' operations. Other operations are all ignored. The length of CIGAR_Str should be always less than 100 bytes or "-1" is returned.
// Staring_Points and Section_Length are empty arrays to write the sections. The minimum length of each array is 6 items.
// The length of a section is its length on the chromosome, namely 'I' is ignored but 'D' is added into the length.
// This function ignores all sections from the 7-th.

// This function returns the number of sections found in the CIGAR string. It returns -1 if the CIGAR string cannot be parsed.

int RSubread_parse_CIGAR_string(const char * CIGAR_Str, unsigned int * Staring_Points, unsigned short * Section_Length);


// This function try to find the attribute value of a given attribute name from the extra column string in GTF/GFF.
// If the value is found, it returns the length of the value (must be > 0 by definition), or -1 if no attribute is found or the format is wrong.

int GTF_extra_column_value(const char * Extra_Col, const char * Target_Name, char * Target_Value, int TargVal_Size);


// Replacing `rep' with `with' in `orig'. 
// Rhe return value must be freed if it is not NULL.
char *str_replace(char *orig, char *rep, char *with) ;


// rule: the string is ABC123XXXXXX...
// // This is the priroity:
// // First, compare the letters part.
// // Second, compare the pure numeric part.
// // Third, compare the remainder.
int strcmp_number(char * s1, char * s2);

#endif
