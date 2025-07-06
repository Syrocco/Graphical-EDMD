CC = gcc
G ?= 0
SANITIZE ?= 0

ifeq ($(G), 0)
	CFLAGS = -Ofast -Wall -Wextra -lm -fopenmp
else
	CFLAGS = -Ofast graphics.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
endif

ifeq ($(SANITIZE), 1)
	CFLAGS += -fsanitize=address -g
endif

EDMD: EDMD.c parser.o
	$(CC) EDMD.c parser.c $(CFLAGS) $(EDMD_FLAGS) -DG=$(G)