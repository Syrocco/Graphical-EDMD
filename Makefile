# Detect macOS and use Homebrew gcc if on Mac
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# Find the Homebrew gcc version (e.g., gcc-14, gcc-15, etc.)
	CC := $(shell ls /opt/homebrew/bin/gcc-* 2>/dev/null | head -1 | xargs basename 2>/dev/null)
	ifeq ($(CC),)
		# Fallback: try /usr/local/bin for Intel Macs
		CC := $(shell ls /usr/local/bin/gcc-* 2>/dev/null | head -1 | xargs basename 2>/dev/null)
	endif
	ifeq ($(CC),)
		# No real GCC found, use clang but warn about potential issues
		CC = clang
		$(warning Warning: No GCC found, using clang. OpenMP support may be limited.)
	endif
else
	CC = gcc
endif


# Check if OpenMP is supported by the compiler
# For clang on macOS, might need libomp installed via: brew install libomp
OPENMP_SUPPORT := $(shell echo 'int main(){return 0;}' | $(CC) -fopenmp -x c - -o /tmp/openmp_test 2>/dev/null && echo "yes" || echo "no"; rm -f /tmp/openmp_test 2>/dev/null)
ifeq ($(OPENMP_SUPPORT),no)
	ifneq ($(findstring clang,$(CC)),)
		# Try clang with explicit libomp path (common on macOS with Homebrew)
		OPENMP_SUPPORT := $(shell echo 'int main(){return 0;}' | $(CC) -Xpreprocessor -fopenmp -lomp -x c - -o /tmp/openmp_test 2>/dev/null && echo "clang" || echo "no"; rm -f /tmp/openmp_test 2>/dev/null)
	endif
endif


# Source files in src directory
SRC_DIR = src
SOURCES = $(SRC_DIR)/EDMD.c $(SRC_DIR)/parser.c $(SRC_DIR)/boop.c $(SRC_DIR)/pcf.c $(SRC_DIR)/voronoi_edmd.c
GRAPHICS_SOURCES = $(SOURCES) $(SRC_DIR)/graphics.c

# Include directory for headers
INCLUDE_FLAGS = -I$(SRC_DIR)

ifeq ($(OPENMP_SUPPORT),yes)
	OPENMP_FLAGS = -fopenmp
else
	OPENMP_FLAGS = 
endif

# Detect OpenMP flags for clang
ifeq ($(OPENMP_SUPPORT),clang)
	OPENMP_FLAGS = -Xpreprocessor -fopenmp -lomp
endif

# Platform-specific graphics libraries
ifeq ($(UNAME_S),Darwin)
	# macOS - use the exact same flags as the working raylib Makefile
	RAYLIB_PREFIX = $(shell brew --prefix raylib)
	GRAPHICS_LIBS = -L$(RAYLIB_PREFIX)/lib -lraylib -framework OpenGL -framework Cocoa -framework IOKit -framework CoreAudio -framework CoreVideo
	GRAPHICS_INCLUDES = -I$(RAYLIB_PREFIX)/include
else
	# Linux - use X11 and GL libraries
	GRAPHICS_LIBS = -lraylib -lGL -lm  -ldl -lrt -lX11
	GRAPHICS_INCLUDES =
endif

EDMD: $(SOURCES)
	$(CC) $(SOURCES) $(INCLUDE_FLAGS) -O3 -ffast-math -Wall -Wextra -lm $(OPENMP_FLAGS) -DG=0

# Graphics target - optimized build with graphics libraries
graphics: $(GRAPHICS_SOURCES)
	$(CC) $(GRAPHICS_SOURCES) $(INCLUDE_FLAGS) $(GRAPHICS_INCLUDES) -std=c99 -Wall -Wno-missing-braces -Wno-unused-value -O3 -ffast-math $(GRAPHICS_LIBS) -lm $(OPENMP_FLAGS) -DG=1

# Debug target - no optimization, with debug symbols
debug: $(SOURCES)
	$(CC) $(SOURCES) $(INCLUDE_FLAGS) -O0 -g -Wall -Wextra -lm -DG=0

# Clean target
clean:
	rm -f EDMD EDMD_debug graphics_build *.o

# Make EDMD the default target
.DEFAULT_GOAL := EDMD