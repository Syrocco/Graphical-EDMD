CC = gcc
G ?= 1


ifeq ($(G), 0)
	CFLAGS = -Ofast -Wall -Wextra -lm 
else
	CFLAGS = -fno-omit-frame-pointer -Ofast graphics.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -fopenmp
endif

EDMD: EDMD.c
	$(CC) EDMD.c $(CFLAGS) -g -DG=$(G)

