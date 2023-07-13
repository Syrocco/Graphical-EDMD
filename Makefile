all: EDMD

EDMD: EDMD.c EDMD.h
	gcc EDMD.c -Ofast -Wall -march=native -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -g
