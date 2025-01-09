/*
 *	Author: Etienne Fayen
 *
*/
#include "parser.h"

/*** Open/close dump file ***/

Dump* dump_open(char* path,char mode){
	/*
	 * Creates a dump object.
	 * Initialises path and mode.
	 * Calls the proper initialisation function depending on the mode.
	 */
	Dump* dump=malloc(sizeof(Dump));
	if(dump==NULL){
		printf("[dump_open]: Failed to allocate dump. Exit.\n");
		exit(1);
	}
	dump->path=path;
	dump->mode=mode;
	if(mode=='r'){init_rdump(dump);}
	else if(mode=='w'){init_wdump(dump);}
	else{
		printf("[dump_open]: Unknown mode '%c'\n",mode);
		exit(1);
	}
	return dump;
}

void init_rdump(Dump* dump){
	/* 
	 * Initialises the input stream dump->file.
	 * Initialises dump attributes to default values.
	 * Calls jump_to_frame(0,dump) to initialise frame specific variables.
	 */
	dump->file=fopen(dump->path,"r");
	if(dump->file == NULL){
		printf("[init_rdump]: Failed to open %s. Exit.\n",dump->path);
		exit(1);
	}
	//Read the file by blocks of size BLOCKSIZE and count the number of TIMESTEP sections
	dump->nframes=0;
	int nframemax=1000;   //Will grow automatically if needed
	int growthfactor=2;
	dump->framepos=(long int*)malloc(nframemax*sizeof(long int));
	if(dump->framepos==NULL){
		printf("[init_rdump]: Failed to allocate an array of size %lu. Exit.\n",nframemax*sizeof(long int));
		exit(1);
	}
	char* block=(char*)malloc(BLOCKSIZE*sizeof(char));
	if(block==NULL){
		printf("[init_rdump]: Failed to allocate a buffer of size %d. Exit.\n",BLOCKSIZE);
		exit(1);
	}
	int nblocks=0;  //Counts how many blocks have been searched to keep tracks of the position of the frame in the stream.
	int overlap=strlen("ITEM: TIMESTEP")-1;  //Blocs overlap by this amount to avoid missing frames where "ITEM: TIMESTEP" is at a bloc boundary
	while(feof(dump->file)==0){
		//Null terminate the block for strstr
		block[fread(block,sizeof(char),BLOCKSIZE-1,dump->file)]='\0';
		char* framestart=block;
		while((framestart=strstr(framestart,"ITEM: TIMESTEP"))!=NULL){
			if(dump->nframes==nframemax){
				nframemax*=growthfactor;
				dump->framepos=(long int *)realloc(dump->framepos,nframemax*sizeof(long int));
				if(dump->framepos==NULL){
					printf("[init_rdump]: Failed to allocate an array of size %lu. Exit.\n",nframemax*sizeof(long int));
					exit(1);
				}
			}
			dump->framepos[dump->nframes]=nblocks*(BLOCKSIZE-1-overlap)+framestart-block;
			++dump->nframes;
			++framestart;   //Offset by one to not find the same frame twice
		}
		if(feof(dump->file)==0){
			fseek(dump->file,-overlap,SEEK_CUR); //Move file cursor to overlap blocks
		}
		++nblocks;
	}
	free(block);
	dump->framepos=(long int*)realloc(dump->framepos,(dump->nframes+1)*sizeof(long int));
	if(dump->framepos==NULL){
		printf("[init_rdump]: Failed to allocate an array of size %lu. Exit.\n",nframemax*sizeof(long int));
		exit(1);
	}
	dump->framepos[dump->nframes]=ftell(dump->file);  //Save eof position
	rewind(dump->file);
	dump->curframe=0;
	dump->framesize=0;
	dump->framebuff=NULL;
	jump_to_frame(0,dump);
}

void init_wdump(Dump* dump){
	/*
	 * Creates a dump object in write mode and initialises the stream dump->file to w+ mode.
	 * Initialises nframes, framepos and curframe.
	 */
	dump->file=fopen(dump->path,"w+");
	if(dump->file == NULL){
		printf("[init_wdump]: Could not open %s.\n",dump->path);
		exit(1);
	}
	dump->sizeframepos=100;   // Initial size of framepos array
	dump->framepos=malloc(dump->sizeframepos*sizeof(long int));
	if(dump->framepos==NULL){
		printf("[init_wdump]: Could not allocate framepos array of size %d.\n",dump->sizeframepos);
		exit(1);
	}
	dump->nframes=0;
	dump->curframe=0;
	dump->framesize=0;
	dump->framebuff=malloc(1*sizeof(char));
	dump->framebuff[0]='\0';
	dump->nsec=0;
}

void dump_close(Dump* dump){
	fclose(dump->file);
	free(dump->framebuff);
	free(dump->framepos);
	free(dump);
}

/*** Frame selection ***/

void next_frame(Dump* dump){
	jump_to_frame(dump->curframe+1,dump);
}

void jump_to_frame(int i,Dump* dump){
	/*
	 * Initialises frame specific variables for the specified frame.
	 */
	if(i>=dump->nframes){
		printf("[jump_to_frame]: frame %d out of bounds (nframes = %d)\n",i,dump->nframes);
		exit(1);
	}
	dump->framesize=dump->framepos[i+1]-dump->framepos[i];
	free(dump->framebuff);
	dump->framebuff=(char*)malloc((dump->framesize+1)*sizeof(char));
	if(dump->framebuff==NULL){
		printf("[jump_to_frame]: Failed to assign frame buffer of size %lu.Exit.\n",(dump->framesize+1)*sizeof(char));
		exit(1);
	}
	fseek(dump->file,dump->framepos[i],SEEK_SET);
	int num = fread(dump->framebuff,sizeof(char),dump->framesize,dump->file);
	if (num == 0){
		printf("fread error");
	}
	dump->framebuff[dump->framesize]='\0';
	dump->curframe=i;
	find_sections(dump);
}

void find_sections(Dump* dump){
	/*
	 * Lists sections present in the current framebuff.
	 * Saves their labels in seclabels and pointer to their beginning in framebuff at the same index in secpos.
	 * Saves the number of sections in nsec.
	 */
	dump->nsec=0;
	char* startsec=dump->framebuff;
	while((startsec=strstr(startsec,"ITEM:"))!=NULL){
		if(dump->nsec==NSECMAX){
			printf("[find_sections]: Maximum number of section reached. Recompile the parser with a larger NSECMAX.\n");
			exit(1);
		}
		dump->secpos[dump->nsec]=startsec;
		int lenlabel=0;
		while(startsec[strlen("ITEM: ")+lenlabel]!=' '&&startsec[strlen("ITEM: ")+lenlabel]!='\n'){++lenlabel;}
		if(lenlabel+1>LINESIZE){
			printf("[find_sections]: Section label is longer than LINESIZE.\n");
			exit(1);
		}
		strlcpy(dump->seclabels[dump->nsec],startsec+strlen("ITEM: "),lenlabel+1);
		++dump->nsec;
		++startsec;
	}
	dump->secpos[dump->nsec]=dump->framebuff+strlen(dump->framebuff);
}


/*** Frame parsing ***/

//Generic functions
int sec_index(char*sec,Dump* dump){
	/*
	 * Returns the index of the section sec in seclabels and secpos for the currently loaded frame.
	 * Returns -1 if the sec is not present in seclabels.
	 * The provided sec is compared against section labels only up to the length of the section labels. This means that "BOXTRUC" will match a section "BOX".
	 */
	for(int i=0;i<dump->nsec;++i){
		if(strncmp(sec,dump->seclabels[i],strlen(dump->seclabels[i]))==0){
			return i;
		}
	}
	return -1;
}

void get_header(char* sec,char* header,int headsize,Dump* dump){
	/*
	 * Puts header (first line) of the section sec of the current loaded frame into the provided header string of size headsize. 
	 * If either the section is missing or the header does not fit into the provided string, displays a message and exits.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_header]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int headlen=0;
	while(dump->secpos[secid][headlen]!='\n'){++headlen;}
	if(headlen+1>headsize){
		printf("[get_header]: The header of section %s in frame %d does not fit in the provided buffer of size %d.\n",dump->seclabels[secid],dump->curframe,headsize);
		exit(1);
	}
	strlcpy(header,dump->secpos[secid],headlen+1);
}

int read_intval(char* sec,int n,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the one with index (n,m) as an integer value.
	 * If (n,m) is out of bounds, prints an error message and exits.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	char secbuff[seclen+1];
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* line=strtok(secbuff,"\n");
	line=strtok(NULL,"\n");   //Skip header line
	for(int i=0;i<n;++i){
		line=strtok(NULL,"\n");
	}
	if(line==NULL){
		printf("[read_intval]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],n);
		exit(1);
	}
	char* tok=strtok(line," \t");
	for(int i=0;i<m;++i){
		tok=strtok(NULL," \t");
	}
	if(tok==NULL){
		printf("[read_intval]: In frame %d, section %s: column %d is out of bounds at line %d.\n",dump->curframe,dump->seclabels[secid],m,n);
		exit(1);
	}
	int val=0;
	sscanf(tok,"%d",&val);
	return val;
}

double read_doubleval(char* sec,int n,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the one with index (n,m) as an integer value.
	 * If (n,m) is out of bounds, prints an error message and exits.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	char secbuff[seclen+1];
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* line=strtok(secbuff,"\n");
	line=strtok(NULL,"\n");   //Skip header line
	for(int i=0;i<n;++i){
		line=strtok(NULL,"\n");
	}
	if(line==NULL){
		printf("[read_intval]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],n);
		exit(1);
	}
	char* tok=strtok(line," \t");
	for(int i=0;i<m;++i){
		tok=strtok(NULL," \t");
	}
	if(tok==NULL){
		printf("[read_intval]: In frame %d, section %s: column %d is out of bounds at line %d.\n",dump->curframe,dump->seclabels[secid],m,n);
		exit(1);
	}
	double val=0;
	sscanf(tok,"%lf",&val);
	return val;
}

float read_floatval(char* sec,int n,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the one with index (n,m) as an integer value.
	 * If (n,m) is out of bounds, prints an error message and exits.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	char secbuff[seclen+1];
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* line=strtok(secbuff,"\n");
	line=strtok(NULL,"\n");   //Skip header line
	for(int i=0;i<n;++i){
		line=strtok(NULL,"\n");
	}
	if(line==NULL){
		printf("[read_intval]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],n);
		exit(1);
	}
	char* tok=strtok(line," \t");
	for(int i=0;i<m;++i){
		tok=strtok(NULL," \t");
	}
	if(tok==NULL){
		printf("[read_intval]: In frame %d, section %s: column %d is out of bounds at line %d.\n",dump->curframe,dump->seclabels[secid],m,n);
		exit(1);
	}
	float val=0;
	sscanf(tok,"%f",&val);
	return val;
}


void read_intcol(char* sec,int* vec,int N,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the column  with index m (shape N*1) as an integer array saved in vec.
	 * If N or m is out of bounds, prints an error message and exits.
	 * vec is assumed to be large enough to contain the N values.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	char* secbuff = calloc(seclen+1, sizeof(char));
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* saveptrline=NULL;
	char* saveptrtok=NULL;
	char* line=strtok_r(secbuff,"\n",&saveptrline);
	for(int i=0;i<N;++i){
		line=strtok_r(NULL,"\n",&saveptrline);
		if(line==NULL){
			printf("[read_intcol]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i);
			exit(1);
		}
		char* tok=strtok_r(line," \t",&saveptrtok);
		for(int j=0;j<m;++j){tok=strtok_r(NULL," \t",&saveptrtok);}
		if(tok==NULL){
			printf("[read_intcol]: In frame %d, section %s, line %d: column %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i,m);
			exit(1);
		}
		sscanf(tok,"%d",&vec[i]);		
	}
	free(secbuff);
}

void read_doublecol(char* sec,double* vec,int N,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the column  with index m (shape N*1) as a double array saved in vec.
	 * If N or m is out of bounds, prints an error message and exits.
	 * vec is assumed to be large enough to contain the N values.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	
	char* secbuff = calloc(seclen+1, sizeof(char));
	
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* saveptrline=NULL;
	char* saveptrtok=NULL;
	char* line=strtok_r(secbuff,"\n",&saveptrline);
	for(int i=0;i<N;++i){
		line=strtok_r(NULL,"\n",&saveptrline);
		if(line==NULL){
			printf("[read_intcol]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i);
			exit(1);
		}
		char* tok=strtok_r(line," \t",&saveptrtok);
		for(int j=0;j<m;++j){tok=strtok_r(NULL," \t",&saveptrtok);}
		if(tok==NULL){
			printf("[read_intcol]: In frame %d, section %s, line %d: column %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i,m);
			exit(1);
		}
		//printf("%d %s\n", i, tok);
		sscanf(tok,"%lf",&vec[i]);		
	}
	free(secbuff);
}

void read_floatcol(char* sec,float* vec,int N,int m,Dump* dump){
	/*
	 * Assuming that section sec contains a N*M table of values, returns the column  with index m (shape N*1) as a float array saved in vec.
	 * If N or m is out of bounds, prints an error message and exits.
	 * vec is assumed to be large enough to contain the N values.
	 */
	int secid=sec_index(sec,dump);
	if(secid<0){
		printf("[get_intval]: No section named %s in frame %d.\n",sec,dump->curframe);
		exit(0);
	}
	int seclen=dump->secpos[secid+1]-dump->secpos[secid];
	
	char* secbuff = calloc(seclen+1, sizeof(char));
	
	strlcpy(secbuff,dump->secpos[secid],seclen+1);
	char* saveptrline=NULL;
	char* saveptrtok=NULL;
	char* line=strtok_r(secbuff,"\n",&saveptrline);
	for(int i=0;i<N;++i){
		line=strtok_r(NULL,"\n",&saveptrline);
		if(line==NULL){
			printf("[read_intcol]: In frame %d, section %s: line %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i);
			exit(1);
		}
		char* tok=strtok_r(line," \t",&saveptrtok);
		for(int j=0;j<m;++j){tok=strtok_r(NULL," \t",&saveptrtok);}
		if(tok==NULL){
			printf("[read_intcol]: In frame %d, section %s, line %d: column %d is out of bounds.\n",dump->curframe,dump->seclabels[secid],i,m);
			exit(1);
		}
		//printf("%d %s\n", i, tok);
		sscanf(tok,"%f",&vec[i]);		
	}
	free(secbuff);
}

//High-level parsing functions

double get_timestep(Dump* dump){
	return read_doubleval("TIMESTEP",0,0,dump);
}

int get_natoms(Dump* dump){
	return read_intval("NUMBER",0,0,dump);
}

double get_boxx(char* h,size_t sh,Dump* dump){
	/*
	 * Returns the box length along the x direction in current frame.
	 * If h is not NULL, puts at most sh bytes of the corresponding header parameter in h.
	 */

	if(h!=NULL){
		char header[LINESIZE];
		get_header("BOX",header,LINESIZE,dump);
		char* tok=strtok(header," \t"); //ITEM:
		tok=strtok(NULL," \t");			//BOX
		tok=strtok(NULL," \t");			//BOUNDS
		tok=strtok(NULL," \t");			//What we want
		strlcpy(h,tok,sh);
	}
	return read_doubleval("BOX",0,1,dump);
}

double get_boxy(char* h,size_t sh,Dump* dump){
	/*
	 * Returns the box length along the y direction in current frame.
	 * If h is not NULL, puts at most sh bytes of the corresponding header parameter in h.
	 */

	if(h!=NULL){
		char header[LINESIZE];
		get_header("BOX",header,LINESIZE,dump);
		char* tok=strtok(header," \t"); //ITEM:
		tok=strtok(NULL," \t");			//BOX
		tok=strtok(NULL," \t");			//BOUNDS
		tok=strtok(NULL," \t");			//xlabel
		tok=strtok(NULL," \t");			//What we want
		strlcpy(h,tok,sh);
	}
	return read_doubleval("BOX",1,1,dump);
}

void get_intatomprop(char* prop,int* ptr,int natoms,Dump* dump){
	/*
	 * Parses the first natoms values of column prop of ITEM: ATOMS as int and stores them in ptr.
	 * Puts the corres
	 * If prop does not exist in atoms property list, prints a message and exit.
	 * natoms is assumed to be smaller or equal than the actual number of atoms (no check performed).
	 * ptr is assumed to be large enough to contain natoms int.
	 */
	int i=atomprop_index(prop,dump);
	if(i<0){
		printf("[get_atompropi]: No atom property %s in frame %d.\n",prop,dump->curframe);
		exit(1);
	}
	read_intcol("ATOMS",ptr,natoms,i,dump);
}

void get_doubleatomprop(char* prop,double* ptr,int natoms,Dump* dump){
	/*
	 * Parses the first natoms values of column prop of ITEM: ATOMS as doubles and stores them in ptr.
	 * If prop does not exist in atoms property list, prints a message and exit.
	 * natoms is assumed to be smaller or equal than the actual number of atoms (no check performed).
	 * ptr is assumed to be large enough to contain natoms int.
	 */
	int i=atomprop_index(prop,dump);
	if(i<0){
		printf("[get_atompropd]: No atom property %s in frame %d.\n",prop,dump->curframe);
		exit(1);
	}
	read_doublecol("ATOMS",ptr,natoms,i,dump);
}

void get_floatatomprop(char* prop,float* ptr,int natoms,Dump* dump){
	/*
	 * Parses the first natoms values of column prop of ITEM: ATOMS as floats and stores them in ptr.
	 * If prop does not exist in atoms property list, prints a message and exit.
	 * natoms is assumed to be smaller or equal than the actual number of atoms (no check performed).
	 * ptr is assumed to be large enough to contain natoms int.
	 */
	int i=atomprop_index(prop,dump);
	if(i<0){
		printf("[get_atompropd]: No atom property %s in frame %d.\n",prop,dump->curframe);
		exit(1);
	}
	read_floatcol("ATOMS",ptr,natoms,i,dump);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"

void get_neigh(int neighmax,int** neigh,Dump* dump){
	/*
	 * Reads section NEIGHBOURS of dump and stores the neighbours lists of atom i in the ith line of neigh, assumed to be large enough to contain all the neighbours lists.
	 * neighmax is the maximum possible number of neighbours for one particle.
	 * A property nb that counts the number of neighbours for each particle must be present in the ATOMS section.
	 */
	int secid=sec_index("NEIGHBOURS",dump);
	if(secid<0){
		printf("[get_neigh]: No section NEIGHBOURS in frame %d.\n",dump->curframe);
		exit(1);
	}
	int natoms=get_natoms(dump);
	int nb[natoms];
	get_intatomprop("nb",nb,natoms,dump);
	int secsize=dump->secpos[secid+1]-dump->secpos[secid];
	char secbuff[secsize+1];
	strlcpy(secbuff,dump->secpos[secid],secsize+1);
	char* ptrline=NULL;
	char* ptrtok=NULL;
	char* line=strtok_r(secbuff,"\n",&ptrline);
	int i=0;
	while((line=strtok_r(NULL,"\n",&ptrline))!=NULL){
		char* tok=strtok_r(line," ",&ptrtok);
		sscanf(tok,"%d:",&i);
		if(nb[i]>neighmax){
			printf("[get_neigh]: The provided table of width %d is too small to hold all %d neighbours id of atom %d.\n",neighmax,nb[i],i);
			exit(1);
		}
		for(int k=0;k<nb[i];++k){
			tok=strtok_r(NULL," ",&ptrtok);
			if(tok==NULL){
				printf("[get_neigh]: ATOMS section states that atom %d has %d neighbours, but only %d could be parsed from the NEIGHBOUR section. Is there a problem in atoms indexing ?\n",i,nb[i],k);
				exit(1);
			}
			sscanf(tok,"%d",&neigh[i*neighmax + k]);
		}
	}
}

#pragma GCC diagnostic pop

// Quiery functions

bool sec_exists(char* sec,Dump* dump){
	/*
	 * Returns 1 if sec exists in seclabels, 0 otherwise.
	 */
	for(int i=0;i<dump->nsec;++i){
		if(strcmp(sec,dump->seclabels[i])==0){return 1;}
	}
	return 0;
}

bool atomprop_exists(char* prop,Dump* dump){
	/*
	 * Returns 1 if prop is found in ATOMS header, 0 otherwise.
	 */
	if(atomprop_index(prop,dump)>=0){return 1;}
	else{return 0;}
}


/*** Writing functions ***/

//Generic functions

void cp_framebuff(Dump* dumpr,Dump* dumpw){
	/*
	 * Copies the current frame buffer of dumpr to dumpw.
	 * Issues a warning message if the buffer is not empty, which indicates that the previous content of the buffer has not been written to the file.
	 */
	if(dumpr->mode!='r'){
		printf("[cp_framebuff]: %s is not opened in read mode.\n",dumpr->path);
		exit(1);
	}
	if(dumpw->mode!='w'){
		printf("[cp_framebuff]: %s is not opened in write mode.\n",dumpw->path);
		exit(1);
	}
	if(strlen(dumpw->framebuff)!=0){
		printf("[cp_framebuff]: Frame buffer for frame %d in write-dump is not empty. Did you forget to write the buffer with write_frame function ? Overwriting anyway.\n",dumpw->nframes);
	}
	dumpw->framebuff=realloc(dumpw->framebuff,(dumpr->framesize+1)*sizeof(char));
	if(dumpw->framebuff==NULL){
		printf("[cp_framebuff]: Failed to allocate frame buffer of size %ld.\n",dumpr->framesize+1);
		exit(1);
	}
	dumpw->framesize=dumpr->framesize;
	strlcpy(dumpw->framebuff,dumpr->framebuff,dumpw->framesize+1);
	find_sections(dumpw);
}

void write_frame(Dump* dumpw){
	/*
	 * Writes the current framebuff to dumpw->file.
	 * Empties the buffer, sets framesize and nsec to 0.
	 * If the buffer is already empty, does nothing.
	 */
	if(dumpw->mode!='w'){
		printf("[write_frame]: dump is not opened in write mode.\n");
		exit(1);
	}
	if(dumpw->framesize==0){return;}
	fseek(dumpw->file,dumpw->framepos[dumpw->nframes],SEEK_SET);
	if((int)fwrite(dumpw->framebuff,sizeof(char),dumpw->framesize,dumpw->file)!=dumpw->framesize){
		printf("[write_frame]: number of written bytes does not correspond to the size of the frame buffer.");
	}
	++dumpw->nframes;
	dumpw->curframe=dumpw->nframes;
	if(dumpw->nframes+1>dumpw->sizeframepos){
		dumpw->framepos=realloc(dumpw->framepos,2*dumpw->sizeframepos*sizeof(long int));
		dumpw->sizeframepos*=2;
		if(dumpw->framepos==NULL){
			printf("[write_frame]: Failed to reallocate framepos array of size %d.\n",dumpw->sizeframepos);
			exit(1);
		}
	}
	dumpw->framepos[dumpw->nframes]=ftell(dumpw->file);
	dumpw->framebuff[0]='\0';
	dumpw->framesize=0;
	dumpw->nsec=0;
}

//High-level writing functions

void cpcurframe(Dump* dumpin,Dump* dumpout){
	/*
	 * Copies the current frame of dumpin to dumpout. 
	 * Can serve as a minimal template for more elaborated writing functions.
	 */
	cp_framebuff(dumpin,dumpout);
	//Any processing of dumpout framebuff can be performed here. Just remember to keep frame size up to date and the buffer null-terminated.
	write_frame(dumpout);
}

void add_intatomprop(int* vals,char* prop,Dump* dumpout,Dump* dumpin){
	/*
	 * Adds integer values vals of property prop in ATOMS section of dumpout's current frame buffer.
	 * If dumpout's current frame is empty, dumpin's current frame is copied into dumpout and the property is added there.
	 * It is assumed that vals contains one int per particle.
	 * Note that this function does not write the modified buffer to allow for multiple buffer modifications before write is performed.
	 */
	if(dumpout->framesize==0){
		cp_framebuff(dumpin,dumpout);
	}
	if(atomprop_index(prop,dumpout)>=0){
		printf("[add_intatomprop]: Property %s already exists in frame %d of %s. Adding it anyway.\n",prop,dumpout->curframe,dumpout->path);
	}
	int natoms=get_natoms(dumpout);
	int atomsecid=sec_index("ATOMS",dumpout);
	int atomsecsize=dumpout->secpos[atomsecid+1]-dumpout->secpos[atomsecid];
	int newatomsecsize=2*atomsecsize;
	char* atomsbuff=(char*)malloc((newatomsecsize+1)*sizeof(char));
	if(atomsbuff==NULL){
		printf("[add_intatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
		exit(1);
	}
	int charcount=0;
	int writtenchar=0;
	// Add prop to ATOMS header
	char* line=strtok(dumpout->secpos[atomsecid],"\n");
	writtenchar=snprintf(atomsbuff,newatomsecsize+1,"%s %s\n",line,prop);
	// Grow newatomsecsize buffer if it is full
	while(writtenchar+charcount>=newatomsecsize){
		newatomsecsize*=2;
		atomsbuff=realloc(atomsbuff,(newatomsecsize+1)*sizeof(char));
		if(atomsbuff==NULL){
			printf("[add_intatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
			exit(1);
		}
		writtenchar=snprintf(atomsbuff+charcount,newatomsecsize+1,"%s %s\n",line,prop);
	}
	charcount+=writtenchar;
	// Add values vals to each line
	for(int i=0;i<natoms;++i){
		line=strtok(NULL,"\n");
		writtenchar=snprintf(atomsbuff+charcount,newatomsecsize-charcount+1,"%s %d\n",line,vals[i]);
		// Grow newatomsecsize buffer if it is full
		while(writtenchar+charcount>=newatomsecsize){
			newatomsecsize*=2;
			atomsbuff=realloc(atomsbuff,(newatomsecsize+1)*sizeof(char));
			if(atomsbuff==NULL){
				printf("[add_intatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
				exit(1);
			}
			writtenchar=snprintf(atomsbuff+charcount,newatomsecsize-charcount+1,"%s %d\n",line,vals[i]);
		}
		charcount+=writtenchar;
	}
	// Write the new frame buffer
	int prevsecsize=dumpout->secpos[atomsecid]-dumpout->secpos[0];
	int nextsecsize=dumpout->secpos[dumpout->nsec]-dumpout->secpos[atomsecid+1];
	dumpout->framesize=prevsecsize+charcount+nextsecsize;
	char* newframebuff=(char*)malloc((dumpout->framesize+1)*sizeof(char));
	strlcpy(newframebuff,dumpout->framebuff,prevsecsize+1);
	strlcpy(newframebuff+prevsecsize,atomsbuff,charcount+1);
	strlcpy(newframebuff+prevsecsize+charcount,dumpout->secpos[atomsecid+1],nextsecsize+1);
	free(atomsbuff);
	free(dumpout->framebuff);
	dumpout->framebuff=newframebuff;
	// Refresh section variables
	find_sections(dumpout);  
}

void add_doubleatomprop(double* vals,char* prop,Dump* dumpout,Dump* dumpin){
	/*
	 * Adds double values vals of property prop in ATOMS section of dumpout's current frame buffer.
	 * If dumpout's current frame is empty, dumpin's current frame is copied into dumpout and the property is added there.
	 * It is assumed that vals contains one double per particle.
	 * Note that this function does not write the modified buffer to allow for multiple buffer modifications before write is performed.
	 */
	if(dumpout->framesize==0){
		cp_framebuff(dumpin,dumpout);
	}
	if(atomprop_index(prop,dumpout)>=0){
		printf("[add_doubleatomprop]: Property %s already exists in frame %d of %s. Adding it anyway.\n",prop,dumpout->curframe,dumpout->path);
	}
	int natoms=get_natoms(dumpout);
	int atomsecid=sec_index("ATOMS",dumpout);
	int atomsecsize=dumpout->secpos[atomsecid+1]-dumpout->secpos[atomsecid];
	int newatomsecsize=2*atomsecsize;
	char* atomsbuff=(char*)malloc((newatomsecsize+1)*sizeof(char));
	if(atomsbuff==NULL){
		printf("[add_doubleatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
		exit(1);
	}
	int charcount=0;
	int writtenchar=0;
	// Add prop to ATOMS header
	char* line=strtok(dumpout->secpos[atomsecid],"\n");
	writtenchar=snprintf(atomsbuff,newatomsecsize+1,"%s %s\n",line,prop);
	// Grow newatomsecsize buffer if it is full
	while(writtenchar+charcount>=newatomsecsize){
		newatomsecsize*=2;
		atomsbuff=realloc(atomsbuff,(newatomsecsize+1)*sizeof(char));
		if(atomsbuff==NULL){
			printf("[add_doubleatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
			exit(1);
		}
		writtenchar=snprintf(atomsbuff+charcount,newatomsecsize+1,"%s %s\n",line,prop);
	}
	charcount+=writtenchar;
	// Add values vals to each line
	for(int i=0;i<natoms;++i){
		line=strtok(NULL,"\n");
		writtenchar=snprintf(atomsbuff+charcount,newatomsecsize-charcount+1,"%s %lf\n",line,vals[i]);
		// Grow newatomsecsize buffer if it is full
		while(writtenchar+charcount>=newatomsecsize){
			newatomsecsize*=2;
			atomsbuff=realloc(atomsbuff,(newatomsecsize+1)*sizeof(char));
			if(atomsbuff==NULL){
				printf("[add_doubleatomprop]: Could not allocate buffer of size %d.\n",newatomsecsize+1);
				exit(1);
			}
			writtenchar=snprintf(atomsbuff+charcount,newatomsecsize-charcount+1,"%s %lf\n",line,vals[i]);
		}
		charcount+=writtenchar;
	}
	// Write the new frame buffer
	int prevsecsize=dumpout->secpos[atomsecid]-dumpout->secpos[0];
	int nextsecsize=dumpout->secpos[dumpout->nsec]-dumpout->secpos[atomsecid+1];
	dumpout->framesize=prevsecsize+charcount+nextsecsize;
	char* newframebuff=(char*)malloc((dumpout->framesize+1)*sizeof(char));
	strlcpy(newframebuff,dumpout->framebuff,prevsecsize+1);
	strlcpy(newframebuff+prevsecsize,atomsbuff,charcount+1);
	strlcpy(newframebuff+prevsecsize+charcount,dumpout->secpos[atomsecid+1],nextsecsize+1);
	free(atomsbuff);
	free(dumpout->framebuff);
	dumpout->framebuff=newframebuff;
	// Refresh section variables
	find_sections(dumpout);  
}

void del_atomprop(char* prop,Dump* dumpout,Dump* dumpin){
	/*
	 * Deletes the column corresponding to property prop in the ATOMS section of dumpout's frame buffer.
	 * If dumpout's frame buffer is empty, copies dumpin's one.
	 * Note that this function does not write the modified buffer to allow for multiple buffer modifications before write is performed.
	 */
	if(dumpout->framesize==0){
		cp_framebuff(dumpin,dumpout);
	}
	int propindex=atomprop_index(prop,dumpout);
	if(propindex<0){
		printf("[del_atomprop]: No property %s in frame %d of %s.\n",prop,dumpout->curframe,dumpout->path);
		exit(1);
	}
	int natoms=get_natoms(dumpout);
	int atomsecid=sec_index("ATOMS",dumpout);
	int atomsecsize=dumpout->secpos[atomsecid+1]-dumpout->secpos[atomsecid];
	char* atomsbuff=(char*)malloc((atomsecsize+1)*sizeof(char));
	if(atomsbuff==NULL){
		printf("[del_atomprop]: Could not allocate buffer of size %d.\n",atomsecsize+1);
		exit(1);
	}
	int writtenchar=0;
	int nprops=propindex+1;
	char* lineptr=NULL;
	char* tokptr=NULL;
	// Remove prop from ATOMS header
	char* line=strtok_r(dumpout->secpos[atomsecid],"\n",&lineptr);
	char* tok=strtok_r(line," ",&tokptr);
	for(int i=0;i<propindex+2;++i){
		writtenchar+=snprintf(atomsbuff+writtenchar,atomsecsize+1-writtenchar,"%s ",tok);
		tok=strtok_r(NULL," ",&tokptr);
	}
	while((tok=strtok_r(NULL," ",&tokptr))!=NULL){
		++nprops;	// Cound properties after prop
		writtenchar+=snprintf(atomsbuff+writtenchar,atomsecsize+1-writtenchar,"%s ",tok);
	}
	atomsbuff[writtenchar-1]='\n';  // Replace last space with new line
	// Remove the corresponding column
	for(int i=0;i<natoms;++i){
		line=strtok_r(NULL,"\n",&lineptr);
		tok=strtok_r(line," ",&tokptr);
		for(int k=0;k<nprops;++k){
			if(k==propindex){tok=strtok_r(NULL," ",&tokptr);continue;}
			else{
				writtenchar+=snprintf(atomsbuff+writtenchar,atomsecsize+1-writtenchar,"%s ",tok);
				tok=strtok_r(NULL," ",&tokptr);
			}
		}
		atomsbuff[writtenchar-1]='\n';  // Replace last space with new line
	}
	// Write the new frame buffer
	int prevsecsize=dumpout->secpos[atomsecid]-dumpout->secpos[0];
	int nextsecsize=dumpout->secpos[dumpout->nsec]-dumpout->secpos[atomsecid+1];
	dumpout->framesize=prevsecsize+writtenchar+nextsecsize;
	char* newframebuff=(char*)malloc((dumpout->framesize+1)*sizeof(char));
	strlcpy(newframebuff,dumpout->framebuff,prevsecsize+1);
	strlcpy(newframebuff+prevsecsize,atomsbuff,writtenchar+1);
	strlcpy(newframebuff+prevsecsize+writtenchar,dumpout->secpos[atomsecid+1],nextsecsize+1);
	free(dumpout->framebuff);
	free(atomsbuff);
	dumpout->framebuff=newframebuff;
	find_sections(dumpout);  // Refresh sections variables
}

void add_section(char* header,int headsize,char* content,int contentsize,Dump* dumpout,Dump* dumpin){
	/*
	 * Adds a section to current framebuff of dumpout, with header of size headsize (excluding null terminating byte), and content of size contentsize (excluding null terminating byte).
	 * If the current buffer of dumpout is empty, copies the one from dumpin before adding the section.
	 * Calls find_sections to update sections variables.
	 */
	if(dumpout->mode!='w'){
		printf("[add_section]: %s is not opened in write mode.\n",dumpout->path);
		exit(1);
	}
	if(dumpin->mode!='r'){
		printf("[add_section]: %s is not opened in read mode.\n",dumpin->path);
		exit(1);
	}
	if(dumpout->framesize==0){
		cp_framebuff(dumpin,dumpout);
	}
	int oldframesize=dumpout->framesize;
	dumpout->framesize+=headsize+contentsize;
	dumpout->framebuff=realloc(dumpout->framebuff,dumpout->framesize+1);
	if(dumpout->framebuff==NULL){
		printf("[add_section]: Failed to allocate frame buffer of size %ld.\n",dumpout->framesize);
		exit(1);
	}
	strlcpy(dumpout->framebuff+oldframesize,header,headsize+1);
	strlcpy(dumpout->framebuff+oldframesize+headsize,content,contentsize+1);
	find_sections(dumpout);
}


/*** Internal functions ***/

int atomprop_index(char* prop,Dump* dump){
	/*
	 * Searches the list of properties in ATOMS header for prop and returns its index.
	 * Returns -1 if prop is not present.
	 * Note that only the len(prop) first char of each property in header is compared against prop. Hence "x" will match with the first encountered property starting with "x".
	 */
	char header[LINESIZE];
	get_header("ATOMS",header,LINESIZE,dump);
	char* tok=strtok(header," \t");   //ITEM:
	tok=strtok(NULL," \t");			  //ATOMS
	int index=0;
	while((tok=strtok(NULL," \t"))!=NULL){
		if(strncmp(tok,prop,strlen(prop))==0){return index;}
		++index;
	}
	return -1;
}


size_t strlcpy(char* dst,const char* src,size_t maxlen){
	const size_t srclen = strlen(src);
	if(srclen+1<maxlen){
		memcpy(dst,src,srclen+1);
	}
	else if(maxlen!=0){
		memcpy(dst,src,maxlen-1);
		dst[maxlen-1]='\0';
	}
	return srclen;
}
