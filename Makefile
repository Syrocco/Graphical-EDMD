CC = gcc
G ?= 1


ifeq ($(G), 0)
	CFLAGS = -Ofast -Wall -Wextra -lm
else
	CFLAGS = -Ofast graphics.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
endif


EDMD: EDMD.c
	$(CC) EDMD.c $(CFLAGS) -DG=$(G)
