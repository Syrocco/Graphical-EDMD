# Detect macOS and use Homebrew gcc if on Mac
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# Find the Homebrew gcc version (e.g., gcc-14, gcc-15, etc.)
	CC := $(shell ls /opt/homebrew/bin/gcc-* 2>/dev/null | head -1 | xargs basename 2>/dev/null || echo gcc)
	ifeq ($(CC),gcc)
		# Fallback: try /usr/local/bin for Intel Macs
		CC := $(shell ls /usr/local/bin/gcc-* 2>/dev/null | head -1 | xargs basename 2>/dev/null || echo gcc)
	endif
else
	CC = gcc
endif

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

# Default target (G=0)
EDMD: EDMD.c parser.o boop.c pcf.c
	$(CC) EDMD.c parser.c boop.c pcf.c $(CFLAGS) $(EDMD_FLAGS) -DG=0

# Graphics target (G=1)
graphics: EDMD.c parser.o boop.c pcf.c
	$(CC) EDMD.c parser.c boop.c pcf.c -Ofast graphics.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 $(EDMD_FLAGS) -DG=1

# Make EDMD the default target
.DEFAULT_GOAL := EDMD