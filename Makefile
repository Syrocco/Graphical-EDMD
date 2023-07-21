CC = gcc
G ?= 1

ifeq ($(G), 0)
	CFLAGS = -Ofast -Wall -Wextra -march=native -lm
else
	CFLAGS = -Ofast -march=native -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -fopenmp
endif

EDMD: EDMD.c
	$(CC) EDMD.c quartic_real.c $(CFLAGS) -g -DG=$(G)

