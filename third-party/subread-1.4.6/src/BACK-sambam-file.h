#ifndef _SAMBAM_FILE_H_
#define _SAMBAM_FILE_H_

typedef unsigned char BS_uint_8;
typedef unsigned short BS_uint_16;
typedef unsigned int BS_uint_32;

typedef struct 
{
	FILE * os_file;
	unsigned long long next_block_header_offset;
	int file_mode;
} BamSam_FILE;

typedef struct
{
	BS_uint_8 ID1;
	BS_uint_8 ID2;
	BS_uint_8 CM;
	BS_uint_8 FLG;
	BS_uint_32 MTIME;
	BS_uint_8 XFL;
	BS_uint_8 OS;
	BS_uint_16 XLEN;
} BGZF_Header;

typedef struct
{
	BS_uint_8 SI1;
	BS_uint_8 SI2;
	BS_uint_16 SLEN;
	BS_uint_16 BSIZE;
} RCF1952_Subfield;

#define SAMBAM_FILE_SAM	10
#define SAMBAM_FILE_BAM 20

// this function returns 0 if OK, or a minus value if the file reaches EOF.
// the file pointer is put to the first byte of the first subfield after the header.
int get_next_BGZF_block_header(BamSam_FILE * fp, BGZF_Header * header);

// The file pointer is put to the first byte of CDATA.
// There should be a BamSam subfield found so the function can return 0; if not, the return value is minus.
int get_RFC1952_subfield(BamSam_FILE * fp, RCF1952_Subfield * field);

// This function moves 8 bytes of fp, putting the file pointer to the first byte of the next BamSam block.
int finalise_BamSam_block(BamSam_FILE * fp);

// This function opens a BamSam file in read-only mode.
// FILE_MODE specifies if it is a SAM file or a BAM file.
// It returns NULL if fopen(fn,"r") == NULL.
BamSam_FILE * fopen_BamSam(const char * fn, int FILE_MODE);

// just like feof(fp)
int feof_BamSam(BamSam_FILE *fp);
#endif
