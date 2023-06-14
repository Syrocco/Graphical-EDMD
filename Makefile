all: EDMD GEN

EDMD: EDMD.c EDMD.h
	gcc EDMD.c -lm -Ofast -Wall
 
GEN: mdgrow.c mdgrow.h
	gcc mdgrow.c -lm -Ofast -o gen -mcmodel=medium
