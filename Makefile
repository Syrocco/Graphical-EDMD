all: EDMD

EDMD: EDMD.c EDMD.h
	gcc EDMD.c -Wall -march=native -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -g -fopenmp -Ofast
