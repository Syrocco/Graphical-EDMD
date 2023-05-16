all: EDMD
EDMD: EDMD.c EDMD.h
	gcc EDMD.c -lm -Ofast -Wall
