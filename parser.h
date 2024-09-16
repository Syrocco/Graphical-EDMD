#ifndef __PARSER_H
#define __PARSER_H

/*
 * A library for parsing LAMMPS-format dump files 
 */
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>

#define BLOCKSIZE 1048576//Size of the buffer for IO access to the file
#define LINESIZE 256     //Max size of line in dump file
#define NSECMAX 10       //Maximum number of sections allowed in one frame

typedef struct Dump{
	/*
	 * Can be opened in read 'r' or write 'w' mode.
	 * In write mode, framepos contains only two values corresponding to the position of the beginning and end of the current frame.
	 * Frame size is 0 in write mode.
	 */
	//File specific variables
	char* path;
	FILE* file;
	char mode;
	//Dump specific variables
	int nframes;    //Number of frames
	long int* framepos;     //Position of the beginning of each frame in the stream. Last value is the position of the end of the file. Array of size nframes+1
	int sizeframepos;		// Size of the framepos array, only used for dump opened in write mode to reallocate a larger array if needed
	//Frame specific variables
	int curframe;	 //Current frame id
	long int framesize;  // Number of chars in frame (excluding terminating null byte)
	char* framebuff; // Buffer containing the current frame (must be able to hold at least framesize+1 char)
	int nsec;		 //Number of sections in the current frame
	char seclabels[NSECMAX][LINESIZE]; //List of the labels of available sections in the current frame
	char* secpos[NSECMAX+1];  //Pointer to the beginning of each section in framebuff. The last pointer corresponds to the end of the frame buffer.
}Dump;

/***** User methods *****/
//Open/close dump file
Dump* dump_open(char* path,char mode);
void init_rdump(Dump* dump);
void init_wdump(Dump* dump);
void dump_close(Dump* dump);

//Frame selection
void next_frame(Dump* dump);
void jump_to_frame(int i,Dump* dump);

//Frame parsing
//Generic functions that can be used to write higher-level parsing functions
int sec_index(char* sec,Dump* dump);
void get_header(char* sec,char* header,int headsize,Dump* dump);
int read_intval(char* sec,int n,int m,Dump* dump);
double read_doubleval(char* sec,int n,int m,Dump* dump);
float read_floatleval(char* sec,int n,int m,Dump* dump);
void read_intcol(char* sec,int* vec,int N,int m,Dump* dump);
void read_doublecol(char* sec,double* vec,int N,int m,Dump* dump);
void read_floatcol(char* sec,float* vec,int N,int m,Dump* dump);
//Higher-level parsing functions
double get_timestep(Dump* dump);
int get_natoms(Dump* dump);
double get_boxx(char* h,size_t sh,Dump* dump);
double get_boxy(char* h,size_t sh,Dump* dump);
void get_intatomprop(char* prop,int* ptr,int natoms,Dump* dump);
void get_doubleatomprop(char* prop,double* ptr,int natoms,Dump* dump);
void get_floatatomprop(char* prop,float* ptr,int natoms,Dump* dump);
void get_neigh(int neighmax, int** neigh,Dump* dump);
//Quiery functions
bool sec_exists(char* sec,Dump* dump);
bool atomprop_exists(char* prop,Dump* dump);

//Writing functions
//Generic functions that can be used to write higher-level writing functions
void cp_framebuff(Dump* dumpr,Dump* dumpw);
void write_frame(Dump* dumpw);
//Higher-level writing functions
void cpcurframe(Dump* dumpin,Dump* dumpout);
void add_intatomprop(int* vals,char* prop,Dump* dumpout,Dump* dumpin);
void add_doubleatomprop(double* vals,char* prop,Dump* dumpout,Dump* dumpin);
void del_atomprop(char* prop,Dump* dumpout,Dump* dumpin);
void add_section(char* header,int headsize,char*content,int contentsize,Dump* dumpout,Dump* dumpin);

//Internal functions, not meant to be called by user
void find_sections(Dump* dump);
int atomprop_index(char* prop,Dump* dump);

//BSD safe copy function
size_t strlcpy(char * dst, const char * src, size_t maxlen);

#endif
