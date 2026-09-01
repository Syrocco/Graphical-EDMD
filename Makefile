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

# The graphical executable contains Dear ImGui/rlImGui C++ sources.  Keep the
# C and C++ compiler families together so that OpenMP and runtime libraries
# agree at the final link step.
ifeq ($(UNAME_S),Darwin)
	ifneq ($(filter gcc-%,$(CC)),)
		CXX := $(patsubst gcc-%,g++-%,$(CC))
	else ifneq ($(findstring clang,$(CC)),)
		CXX := clang++
	else
		CXX := c++
	endif
else
	CXX = g++
endif


OPENMP_SUPPORT := $(shell echo 'int main(){return 0;}' | $(CC) -fopenmp -x c - -o /dev/null 2>/dev/null && echo "yes" || echo "no")
ifeq ($(OPENMP_SUPPORT),no)
	ifneq ($(findstring clang,$(CC)),)
		# Try clang with explicit libomp path (common on macOS with Homebrew)
		OPENMP_SUPPORT := $(shell echo 'int main(){return 0;}' | $(CC) -Xpreprocessor -fopenmp -lomp -x c - -o /dev/null 2>/dev/null && echo "clang" || echo "no")
	endif
endif



SRC_DIR = src
SOURCES = $(SRC_DIR)/EDMD.c $(SRC_DIR)/parser.c $(SRC_DIR)/boop.c $(SRC_DIR)/pcf.c $(SRC_DIR)/voronoi_edmd.c $(SRC_DIR)/osmosis.c $(SRC_DIR)/intruder.c $(SRC_DIR)/field.c
GRAPHICS_SOURCES = $(SOURCES) $(SRC_DIR)/graphics.c \
	$(SRC_DIR)/sim_config.c $(SRC_DIR)/sim_lifecycle.c \
	$(SRC_DIR)/live_analysis.c

IMGUI_DIR = $(SRC_DIR)/vendor/imgui
RLIMGUI_DIR = $(SRC_DIR)/vendor/rlimgui
IMGUI_SOURCES = $(IMGUI_DIR)/imgui.cpp \
	$(IMGUI_DIR)/imgui_draw.cpp \
	$(IMGUI_DIR)/imgui_tables.cpp \
	$(IMGUI_DIR)/imgui_widgets.cpp
RLIMGUI_SOURCES = $(RLIMGUI_DIR)/rlImGui.cpp
GRAPHICS_VENDOR_CXX_SOURCES = $(IMGUI_SOURCES) $(RLIMGUI_SOURCES)
# Add project-owned C++ GUI modules here.  Unlike vendored sources, these are
# compiled with the EDMD dimension and graphics feature macros.
GRAPHICS_APPLICATION_CXX_SOURCES = $(SRC_DIR)/gui_imgui.cpp \
	$(SRC_DIR)/gui_panels.cpp $(SRC_DIR)/gui_analysis.cpp

# Select the simulation dimension at build time, e.g. `make graphics THREE_D=1`.
# When omitted, EDMD.h supplies the default so source-level configurations keep
# working as before.
THREE_D ?=
ifneq ($(strip $(THREE_D)),)
DIMENSION_FLAGS = -DTHREE_D=$(THREE_D)
endif

# Compile-time EDMD capabilities.  Keep these definitions identical across
# every project-owned C and C++ translation unit so that configuration/UI code
# cannot disagree with the simulation core about the active feature set.
# Vendored ImGui/rlImGui sources deliberately do not receive these macros.
STRESS ?= 0
HEAT ?= 0
TANGENTIAL ?= 0
INTRUDER ?= 0
EDMD_VERIFY_COLLISION_ACCOUNTING ?= 0
EDMD_VERIFY_QUEUE ?= 0
EDMD_VERIFY_OVERLAPS ?= 0
NONADDITIVE ?= 1
EDMD_FEATURE_FLAGS = -DSTRESS=$(STRESS) -DHEAT=$(HEAT) \
	-DTANGENTIAL=$(TANGENTIAL) -DINTRUDER=$(INTRUDER) \
	-DNONADDITIVE=$(NONADDITIVE) \
	-DEDMD_VERIFY_COLLISION_ACCOUNTING=$(EDMD_VERIFY_COLLISION_ACCOUNTING) \
	-DEDMD_VERIFY_QUEUE=$(EDMD_VERIFY_QUEUE) \
	-DEDMD_VERIFY_OVERLAPS=$(EDMD_VERIFY_OVERLAPS)

# Reject combinations whose source models are intrinsically incompatible.
# Fixed-axis tangential *delta* remains available in 3D with TANGENTIAL=0;
# TANGENTIAL=1 selects the older scalar-spin rough-disk model.
ifeq ($(strip $(THREE_D)),1)
ifneq ($(strip $(TANGENTIAL)),0)
$(error TANGENTIAL=1 is the 2D scalar-spin rough-particle model; use the 3D tangential-delta GUI option instead)
endif
ifneq ($(strip $(INTRUDER)),0)
$(error INTRUDER=1 uses a polygonal 2D collision model and cannot be built with THREE_D=1)
endif
endif
ifneq ($(strip $(TANGENTIAL)),0)
ifneq ($(strip $(STRESS)$(HEAT)),00)
$(error TANGENTIAL=1 does not yet have complete STRESS/HEAT transport accounting)
endif
endif
ifneq ($(strip $(INTRUDER)),0)
ifneq ($(strip $(STRESS)$(HEAT)),00)
$(error INTRUDER=1 does not yet have complete STRESS/HEAT transport accounting)
endif
endif

# Compiler flags are not part of Make's dependency graph.  Isolating each
# dimension and feature combination prevents incompatible objects from being
# reused in a subsequent graphics build.
GRAPHICS_DIMENSION = $(if $(strip $(THREE_D)),$(THREE_D),default)
GRAPHICS_FEATURE_SET = stress$(STRESS)-heat$(HEAT)-tangent$(TANGENTIAL)-intruder$(INTRUDER)-contacts$(NONADDITIVE)-verify$(EDMD_VERIFY_COLLISION_ACCOUNTING)-queue$(EDMD_VERIFY_QUEUE)-overlap$(EDMD_VERIFY_OVERLAPS)
GRAPHICS_BUILD_ROOT = build/graphics
GRAPHICS_BUILD_DIR = $(GRAPHICS_BUILD_ROOT)/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
GRAPHICS_C_OBJECTS = $(patsubst %.c,$(GRAPHICS_BUILD_DIR)/%.o,$(GRAPHICS_SOURCES))
GRAPHICS_APPLICATION_CXX_OBJECTS = $(patsubst %.cpp,$(GRAPHICS_BUILD_DIR)/%.o,$(GRAPHICS_APPLICATION_CXX_SOURCES))
# ImGui and rlImGui are dimension-independent and can be reused by 2D and 3D
# links.  They must not see EDMD's short `G` macro, which collides with an
# internal ImGui color variable.
GRAPHICS_VENDOR_OBJECTS = $(patsubst $(SRC_DIR)/vendor/%.cpp,$(GRAPHICS_BUILD_ROOT)/vendor/%.o,$(GRAPHICS_VENDOR_CXX_SOURCES))
GRAPHICS_OBJECTS = $(GRAPHICS_C_OBJECTS) $(GRAPHICS_APPLICATION_CXX_OBJECTS) $(GRAPHICS_VENDOR_OBJECTS)
GRAPHICS_DEPENDENCIES = $(GRAPHICS_OBJECTS:.o=.d)

# Core-test artifacts use the same build-configuration key as project graphics
# objects. Separate paths prevent concurrent 2D/3D matrix jobs from racing on
# an ABI-incompatible EDMD object or executable.
PREDICTION_TEST_BUILD_DIR = build/tests/prediction/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
PREDICTION_TEST_CORE_OBJECT = $(PREDICTION_TEST_BUILD_DIR)/EDMD_prediction_test.o
PREDICTION_TEST_BINARY = $(PREDICTION_TEST_BUILD_DIR)/test_prediction
CALENDAR_TEST_BUILD_DIR = build/tests/calendar/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
CALENDAR_TEST_CORE_OBJECT = $(CALENDAR_TEST_BUILD_DIR)/EDMD_calendar_test.o
CALENDAR_TEST_BINARY = $(CALENDAR_TEST_BUILD_DIR)/test_event_calendar
LIVE_ANALYSIS_TEST_BUILD_DIR = build/tests/live-analysis/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
LIVE_ANALYSIS_TEST_BINARY = $(LIVE_ANALYSIS_TEST_BUILD_DIR)/test_live_analysis
WALL_PRESSURE_TEST_BUILD_DIR = build/tests/wall-pressure/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
WALL_PRESSURE_TEST_CORE_OBJECT = $(WALL_PRESSURE_TEST_BUILD_DIR)/EDMD_wall_pressure_test.o
WALL_PRESSURE_TEST_BINARY = $(WALL_PRESSURE_TEST_BUILD_DIR)/test_wall_pressure
SPHERICAL_WALL_TEST_BUILD_DIR = build/tests/spherical-wall/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
SPHERICAL_WALL_TEST_CORE_OBJECT = $(SPHERICAL_WALL_TEST_BUILD_DIR)/EDMD_spherical_wall_test.o
SPHERICAL_WALL_TEST_BINARY = $(SPHERICAL_WALL_TEST_BUILD_DIR)/test_spherical_wall
SPHERICAL_FIELD_TEST_BUILD_DIR = build/tests/spherical-field/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
SPHERICAL_FIELD_TEST_CORE_OBJECT = $(SPHERICAL_FIELD_TEST_BUILD_DIR)/EDMD_spherical_field_test.o
SPHERICAL_FIELD_TEST_BINARY = $(SPHERICAL_FIELD_TEST_BUILD_DIR)/test_spherical_field
CONFIG_TEST_BUILD_DIR = build/tests/config/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
CONFIG_TEST_CORE_OBJECT = $(CONFIG_TEST_BUILD_DIR)/EDMD_config_test.o
CONFIG_TEST_BINARY = $(CONFIG_TEST_BUILD_DIR)/test_sim_config
HEX_LATTICE_TEST_BUILD_DIR = build/tests/hex-lattice/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
HEX_LATTICE_TEST_CORE_OBJECT = $(HEX_LATTICE_TEST_BUILD_DIR)/EDMD_hex_lattice_test.o
HEX_LATTICE_TEST_BINARY = $(HEX_LATTICE_TEST_BUILD_DIR)/test_hexagonal_lattice
FULL_BOX_3D_TEST_BUILD_DIR = build/tests/full-box/3d-stress-heat
FULL_BOX_3D_TEST_CORE_OBJECT = $(FULL_BOX_3D_TEST_BUILD_DIR)/EDMD_full_box_test.o
FULL_BOX_3D_TEST_BINARY = $(FULL_BOX_3D_TEST_BUILD_DIR)/test_planar_wall_3d
LEES_EDWARDS_3D_TEST_BUILD_DIR = build/tests/lees-edwards/3d-stress-heat
LEES_EDWARDS_3D_TEST_CORE_OBJECT = $(LEES_EDWARDS_3D_TEST_BUILD_DIR)/EDMD_lees_edwards_test.o
LEES_EDWARDS_3D_TEST_BINARY = $(LEES_EDWARDS_3D_TEST_BUILD_DIR)/test_lees_edwards_3d
RECONFIGURE_TEST_BUILD_DIR = build/tests/reconfigure/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
RECONFIGURE_TEST_CORE_OBJECT = $(RECONFIGURE_TEST_BUILD_DIR)/EDMD_reconfigure_test.o
RECONFIGURE_TEST_BINARY = $(RECONFIGURE_TEST_BUILD_DIR)/test_reconfigure

PYTHON ?= python3

INCLUDE_FLAGS = -I$(SRC_DIR)
# `-ffast-math` otherwise enables finite-math-only assumptions which can fold
# away the validator's and scheduler's explicit NaN/Inf checks.  Retain the
# other fast-math transforms while making those safety checks observable.
FINITE_MATH_GUARD_FLAGS = -fno-finite-math-only

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

GRAPHICS_CPPFLAGS = $(INCLUDE_FLAGS) $(GRAPHICS_INCLUDES) \
	-I$(IMGUI_DIR) -I$(RLIMGUI_DIR) $(DIMENSION_FLAGS) \
	$(EDMD_FEATURE_FLAGS) -DG=1
GRAPHICS_VENDOR_CPPFLAGS = $(INCLUDE_FLAGS) $(GRAPHICS_INCLUDES) \
	-I$(IMGUI_DIR) -I$(RLIMGUI_DIR)
GRAPHICS_CFLAGS = -std=c99 -Wall -Wno-missing-braces \
	-Wno-unused-value -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -pthread
# Font Awesome is optional in rlImGui and its generated headers are not part
# of the deliberately minimal vendored dependency set.
GRAPHICS_CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -pthread \
	-DNO_FONT_AWESOME
GRAPHICS_VENDOR_CXXFLAGS = -std=c++17 -O2 -pthread \
	-DNO_FONT_AWESOME

# Browser build. Keep the Emscripten toolchain and all of its products separate
# from native objects: WebAssembly objects cannot be reused by the desktop
# linker, and the dimension/feature macros are part of the project ABI.
WEB_CC ?= emcc
WEB_CXX ?= em++
WEB_AR ?= emar
WEB_BUILD_ROOT = build/web
WEB_OBJECT_ROOT = $(WEB_BUILD_ROOT)/objects
WEB_BUILD_DIR = $(WEB_OBJECT_ROOT)/$(GRAPHICS_DIMENSION)/$(GRAPHICS_FEATURE_SET)
WEB_EM_CACHE ?= $(abspath $(WEB_BUILD_ROOT)/emscripten-cache)
WEB_RAYLIB_VERSION = 5.5
WEB_RAYLIB_SOURCE ?= $(WEB_BUILD_ROOT)/raylib-$(WEB_RAYLIB_VERSION)
WEB_RAYLIB_LIBRARY = $(WEB_RAYLIB_SOURCE)/src/libraylib.a
WEB_RAYLIB_CONFIG_STAMP = $(WEB_RAYLIB_SOURCE)/src/.edmd-webgl2-release-noaudio
WEB_SHELL ?= web/shell.html
WEB_LAUNCHER ?= web/launcher.html
WEB_OUTPUT_DIR ?= docs
WEB_OUTPUT_BASENAME ?= EDMD
WEB_INTERMEDIATE_HTML = $(WEB_OUTPUT_DIR)/$(WEB_OUTPUT_BASENAME).html
WEB_OUTPUT_HTML = $(WEB_OUTPUT_DIR)/index.html

WEB_C_OBJECTS = $(patsubst %.c,$(WEB_BUILD_DIR)/%.o,$(GRAPHICS_SOURCES))
WEB_APPLICATION_CXX_OBJECTS = $(patsubst %.cpp,$(WEB_BUILD_DIR)/%.o,$(GRAPHICS_APPLICATION_CXX_SOURCES))
WEB_VENDOR_OBJECTS = $(patsubst $(SRC_DIR)/vendor/%.cpp,$(WEB_BUILD_DIR)/vendor/%.o,$(GRAPHICS_VENDOR_CXX_SOURCES))
WEB_OBJECTS = $(WEB_C_OBJECTS) $(WEB_APPLICATION_CXX_OBJECTS) $(WEB_VENDOR_OBJECTS)
WEB_DEPENDENCIES = $(WEB_OBJECTS:.o=.d)

WEB_PLATFORM_CPPFLAGS = $(INCLUDE_FLAGS) -I$(WEB_RAYLIB_SOURCE)/src \
	-I$(IMGUI_DIR) -I$(RLIMGUI_DIR) -D_DEFAULT_SOURCE \
	-DPLATFORM_WEB -DGRAPHICS_API_OPENGL_ES3
WEB_PROJECT_CPPFLAGS = $(WEB_PLATFORM_CPPFLAGS) $(DIMENSION_FLAGS) \
	$(EDMD_FEATURE_FLAGS) -DG=1
# Do not give Dear ImGui/rlImGui EDMD feature macros or `G`: ImGui uses `G`
# as an internal identifier, and the project macro would corrupt that code.
WEB_VENDOR_CPPFLAGS = $(WEB_PLATFORM_CPPFLAGS)
WEB_CFLAGS = -std=gnu99 -Wall -Wno-missing-braces -Wno-unused-value \
	-O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS)
WEB_CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -DNO_FONT_AWESOME
WEB_VENDOR_CXXFLAGS = -std=c++17 -O2 -DNO_FONT_AWESOME
WEB_LINK_FLAGS = -O3 -sUSE_GLFW=3 -sGL_ENABLE_GET_PROC_ADDRESS=1 \
	-sMIN_WEBGL_VERSION=2 -sMAX_WEBGL_VERSION=2 \
	-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 \
	-sSTACK_SIZE=1048576 -sFORCE_FILESYSTEM=1 -sENVIRONMENT=web \
	-sEXPORTED_RUNTIME_METHODS=ccall --shell-file $(WEB_SHELL)

EDMD: $(SOURCES)
	$(CC) $(SOURCES) $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0


graphics: $(GRAPHICS_OBJECTS)
	$(CXX) $(GRAPHICS_OBJECTS) $(GRAPHICS_LIBS) -lm -pthread $(OPENMP_FLAGS) -o a.out

# Validate or fetch the pinned raylib source before compiling anything that
# includes raylib.h. This phony order-only prerequisite also makes parallel
# first builds safe without forcing already-built objects to be recompiled.
.PHONY: web-raylib-source
web-raylib-source:
	./scripts/fetch_raylib_web.sh "$(WEB_RAYLIB_SOURCE)"

$(WEB_RAYLIB_CONFIG_STAMP): | web-raylib-source
	@touch "$@"

# raylib's upstream object names do not encode graphics/audio flags. Clean it
# once when establishing this configuration so a pre-existing ES2 archive can
# never be mistaken for the required WebGL 2/ES3 build.
$(WEB_RAYLIB_LIBRARY): $(WEB_RAYLIB_CONFIG_STAMP)
	@mkdir -p "$(WEB_EM_CACHE)"
	EM_CACHE="$(WEB_EM_CACHE)" emmake $(MAKE) -C "$(WEB_RAYLIB_SOURCE)/src" clean \
		PLATFORM=PLATFORM_WEB GRAPHICS=GRAPHICS_API_OPENGL_ES3 \
		RAYLIB_BUILD_MODE=RELEASE RAYLIB_MODULE_AUDIO=FALSE \
		CC="$(WEB_CC)" AR="$(WEB_AR)"
	EM_CACHE="$(WEB_EM_CACHE)" emmake $(MAKE) -C "$(WEB_RAYLIB_SOURCE)/src" \
		PLATFORM=PLATFORM_WEB GRAPHICS=GRAPHICS_API_OPENGL_ES3 \
		RAYLIB_BUILD_MODE=RELEASE RAYLIB_MODULE_AUDIO=FALSE \
		CC="$(WEB_CC)" AR="$(WEB_AR)"

$(WEB_BUILD_DIR)/%.o: %.c | web-raylib-source
	@mkdir -p $(dir $@) "$(WEB_EM_CACHE)"
	EM_CACHE="$(WEB_EM_CACHE)" $(WEB_CC) $(WEB_PROJECT_CPPFLAGS) \
		$(WEB_CFLAGS) -MMD -MP -c $< -o $@

$(WEB_BUILD_DIR)/%.o: %.cpp | web-raylib-source
	@mkdir -p $(dir $@) "$(WEB_EM_CACHE)"
	EM_CACHE="$(WEB_EM_CACHE)" $(WEB_CXX) $(WEB_PROJECT_CPPFLAGS) \
		$(WEB_CXXFLAGS) -MMD -MP -c $< -o $@

$(WEB_BUILD_DIR)/vendor/%.o: $(SRC_DIR)/vendor/%.cpp | web-raylib-source
	@mkdir -p $(dir $@) "$(WEB_EM_CACHE)"
	EM_CACHE="$(WEB_EM_CACHE)" $(WEB_CXX) $(WEB_VENDOR_CPPFLAGS) \
		$(WEB_VENDOR_CXXFLAGS) -MMD -MP -c $< -o $@

.PHONY: web web-all web-2d web-3d web-variant web-launcher \
	web-raylib web-serve web-clean
web-raylib: $(WEB_RAYLIB_LIBRARY)

# `make web` produces one static-host-ready site containing both dimensions.
# Build the variants sequentially: they share raylib's WebAssembly archive and
# Emscripten cache, while their project objects remain dimension-specific.
web: web-all

web-all: web-raylib
	$(MAKE) web-2d WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"
	$(MAKE) web-3d WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"
	$(MAKE) web-launcher WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"

# Keep independently buildable variants for development and CI. Their output
# layout is identical to the corresponding directory in the combined site.
web-2d: web-raylib
	$(MAKE) web-variant THREE_D=0 \
		WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)/2d" \
		WEB_OUTPUT_BASENAME="$(WEB_OUTPUT_BASENAME)"

web-3d: web-raylib
	$(MAKE) web-variant THREE_D=1 \
		WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)/3d" \
		WEB_OUTPUT_BASENAME="$(WEB_OUTPUT_BASENAME)"

# Linking is intentionally repeated for this phony target so a missing JS or
# Wasm sidecar is always regenerated, even if index.html is still present.
# This low-level target is also useful for a custom, single-variant output:
#   make web-variant THREE_D=1 WEB_OUTPUT_DIR=/tmp/edmd-web
web-variant: $(WEB_RAYLIB_LIBRARY) $(WEB_OBJECTS) $(WEB_SHELL)
	@mkdir -p "$(WEB_OUTPUT_DIR)" "$(WEB_EM_CACHE)"
	EM_CACHE="$(WEB_EM_CACHE)" $(WEB_CXX) $(WEB_OBJECTS) \
		$(WEB_RAYLIB_LIBRARY) $(WEB_LINK_FLAGS) -lm \
		-o "$(WEB_INTERMEDIATE_HTML)"
	mv -f "$(WEB_INTERMEDIATE_HTML)" "$(WEB_OUTPUT_HTML)"

# The root page contains no dimension-specific assets. It redirects to the
# selected app with project-relative URLs and gives file:// users a useful
# explanation instead of leaving them with an opaque Wasm loading failure.
web-launcher: $(WEB_LAUNCHER)
	@mkdir -p "$(WEB_OUTPUT_DIR)"
	cp "$(WEB_LAUNCHER)" "$(WEB_OUTPUT_HTML)"
	@touch "$(WEB_OUTPUT_DIR)/.nojekyll"
	@rm -f -- "$(WEB_OUTPUT_DIR)/EDMD.html" \
		"$(WEB_OUTPUT_DIR)/EDMD.js" "$(WEB_OUTPUT_DIR)/EDMD.wasm" \
		"$(WEB_OUTPUT_DIR)/EDMD.data" "$(WEB_OUTPUT_DIR)/EDMD.js.map" \
		"$(WEB_OUTPUT_DIR)/EDMD.wasm.map"

web-serve: web
	$(PYTHON) scripts/serve_web.py "$(WEB_OUTPUT_DIR)"

# Keep the verified raylib 5.5 source checkout so a clean rebuild does not
# require another download. Only known generated paths are removed.
web-clean:
	@if [ -f "build/web/raylib-5.5/src/Makefile" ]; then \
		$(MAKE) -C "build/web/raylib-5.5/src" clean \
			PLATFORM=PLATFORM_WEB GRAPHICS=GRAPHICS_API_OPENGL_ES3; \
	fi
	rm -rf -- build/web/objects build/web/emscripten-cache
	rm -rf -- docs/2d docs/3d
	rm -f -- docs/.nojekyll docs/index.html docs/EDMD.html docs/EDMD.js \
		docs/EDMD.wasm docs/EDMD.data docs/EDMD.js.map docs/EDMD.wasm.map

$(GRAPHICS_BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(GRAPHICS_CPPFLAGS) $(GRAPHICS_CFLAGS) $(OPENMP_FLAGS) -MMD -MP -c $< -o $@

$(GRAPHICS_BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(GRAPHICS_CPPFLAGS) $(GRAPHICS_CXXFLAGS) -MMD -MP -c $< -o $@

$(GRAPHICS_BUILD_ROOT)/vendor/%.o: $(SRC_DIR)/vendor/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(GRAPHICS_VENDOR_CPPFLAGS) $(GRAPHICS_VENDOR_CXXFLAGS) -MMD -MP -c $< -o $@


debug: $(SOURCES)
	$(CC) $(SOURCES) $(INCLUDE_FLAGS) -O0 -g $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0

# Profiling build (good symbols and call stacks for time profiling)
profile: $(SOURCES)
	$(CC) $(SOURCES) $(INCLUDE_FLAGS) -O1 -g -fno-omit-frame-pointer -fno-inline $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -o EDMD_profile

.PHONY: graphics prediction-test calendar-test live-analysis-test wall-pressure-test spherical-wall-test spherical-field-test config-test hex-lattice-test full-box-3d-test lees-edwards-3d-test reconfigure-test reconfigure-matrix-test graphics-performance-test
prediction-test:
	@mkdir -p $(PREDICTION_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -Dmain=edmd_program_main -c -o $(PREDICTION_TEST_CORE_OBJECT)
	$(CC) tests/test_prediction.c $(PREDICTION_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O3 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -o $(PREDICTION_TEST_BINARY)
	$(PREDICTION_TEST_BINARY)

calendar-test:
	@mkdir -p $(CALENDAR_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -Dmain=edmd_program_main -c -o $(CALENDAR_TEST_CORE_OBJECT)
	$(CC) tests/test_event_calendar.c $(CALENDAR_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O3 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -o $(CALENDAR_TEST_BINARY)
	$(CALENDAR_TEST_BINARY)

live-analysis-test:
	@mkdir -p $(LIVE_ANALYSIS_TEST_BUILD_DIR)
	$(CC) tests/test_live_analysis.c $(SRC_DIR)/live_analysis.c $(INCLUDE_FLAGS) \
		-O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=1 \
		-o $(LIVE_ANALYSIS_TEST_BINARY)
	$(LIVE_ANALYSIS_TEST_BINARY)

wall-pressure-test:
	@mkdir -p $(WALL_PRESSURE_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_WALL_X=1 -DG=0 -Dmain=edmd_program_main -c -o $(WALL_PRESSURE_TEST_CORE_OBJECT)
	$(CC) tests/test_wall_pressure.c $(WALL_PRESSURE_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O3 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_WALL_X=1 -DG=0 -o $(WALL_PRESSURE_TEST_BINARY)
	$(WALL_PRESSURE_TEST_BINARY)

spherical-wall-test:
	@mkdir -p $(SPHERICAL_WALL_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_CIRCULAR_WALL=1 -DG=0 -Dmain=edmd_program_main -c -o $(SPHERICAL_WALL_TEST_CORE_OBJECT)
	$(CC) tests/test_spherical_wall.c $(SPHERICAL_WALL_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O3 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_CIRCULAR_WALL=1 -DG=0 -o $(SPHERICAL_WALL_TEST_BINARY)
	$(SPHERICAL_WALL_TEST_BINARY)

spherical-field-test:
	@mkdir -p $(SPHERICAL_FIELD_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_CIRCULAR_WALL=1 -DEDMD_ADD_FIELD=1 -DG=0 -Dmain=edmd_program_main -c -o $(SPHERICAL_FIELD_TEST_CORE_OBJECT)
	$(CC) tests/test_spherical_wall.c $(SPHERICAL_FIELD_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O3 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_ADD_CIRCULAR_WALL=1 -DEDMD_ADD_FIELD=1 -DG=0 -o $(SPHERICAL_FIELD_TEST_BINARY)
	$(SPHERICAL_FIELD_TEST_BINARY)

config-test:
	@mkdir -p $(CONFIG_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -Dmain=edmd_program_main -c -o $(CONFIG_TEST_CORE_OBJECT)
	$(CC) tests/test_sim_config.c $(CONFIG_TEST_CORE_OBJECT) $(SRC_DIR)/sim_config.c $(SRC_DIR)/sim_lifecycle.c $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -o $(CONFIG_TEST_BINARY)
	$(CONFIG_TEST_BINARY)

hex-lattice-test:
	@mkdir -p $(HEX_LATTICE_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -Dmain=edmd_program_main -c -o $(HEX_LATTICE_TEST_CORE_OBJECT)
	$(CC) tests/test_hexagonal_lattice.c $(HEX_LATTICE_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DG=0 -o $(HEX_LATTICE_TEST_BINARY)
	$(HEX_LATTICE_TEST_BINARY)

full-box-3d-test:
	@mkdir -p $(FULL_BOX_3D_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) -DTHREE_D=1 -DSTRESS=1 -DHEAT=1 -DTANGENTIAL=0 -DINTRUDER=0 -DNONADDITIVE=1 -DEDMD_ADD_WALL_Z=1 -DEDMD_TEST_PREPARED_STATE=1 -DG=0 -Dmain=edmd_program_main -c -o $(FULL_BOX_3D_TEST_CORE_OBJECT)
	$(CC) tests/test_planar_wall_3d.c $(FULL_BOX_3D_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) -DTHREE_D=1 -DSTRESS=1 -DHEAT=1 -DTANGENTIAL=0 -DINTRUDER=0 -DNONADDITIVE=1 -DEDMD_ADD_WALL_Z=1 -DEDMD_TEST_PREPARED_STATE=1 -DG=0 -o $(FULL_BOX_3D_TEST_BINARY)
	$(FULL_BOX_3D_TEST_BINARY)

lees-edwards-3d-test:
	@mkdir -p $(LEES_EDWARDS_3D_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) -O3 -ffast-math $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) -DTHREE_D=1 -DSTRESS=1 -DHEAT=1 -DTANGENTIAL=0 -DINTRUDER=0 -DNONADDITIVE=1 -DEDMD_ADD_SHEAR=1 -DEDMD_TEST_PREPARED_STATE=1 -DG=0 -Dmain=edmd_program_main -c -o $(LEES_EDWARDS_3D_TEST_CORE_OBJECT)
	$(CC) tests/test_lees_edwards_3d.c $(LEES_EDWARDS_3D_TEST_CORE_OBJECT) $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm $(OPENMP_FLAGS) -DTHREE_D=1 -DSTRESS=1 -DHEAT=1 -DTANGENTIAL=0 -DINTRUDER=0 -DNONADDITIVE=1 -DEDMD_ADD_SHEAR=1 -DEDMD_TEST_PREPARED_STATE=1 -DG=0 -o $(LEES_EDWARDS_3D_TEST_BINARY)
	$(LEES_EDWARDS_3D_TEST_BINARY)

reconfigure-test:
	@mkdir -p $(RECONFIGURE_TEST_BUILD_DIR)
	$(CC) $(SRC_DIR)/EDMD.c $(INCLUDE_FLAGS) $(GRAPHICS_INCLUDES) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_TEST_PREPARED_STATE=1 -DG=1 -Dmain=edmd_program_main -c -o $(RECONFIGURE_TEST_CORE_OBJECT)
	$(CC) tests/test_reconfigure.c $(RECONFIGURE_TEST_CORE_OBJECT) $(SRC_DIR)/sim_config.c $(SRC_DIR)/sim_lifecycle.c $(filter-out $(SRC_DIR)/EDMD.c,$(SOURCES)) $(INCLUDE_FLAGS) $(GRAPHICS_INCLUDES) -O2 $(FINITE_MATH_GUARD_FLAGS) -Wall -Wextra -lm -pthread $(OPENMP_FLAGS) $(DIMENSION_FLAGS) $(EDMD_FEATURE_FLAGS) -DEDMD_TEST_PREPARED_STATE=1 -DG=1 -o $(RECONFIGURE_TEST_BINARY)
	$(RECONFIGURE_TEST_BINARY)

reconfigure-matrix-test:
	$(MAKE) reconfigure-test THREE_D=0 STRESS=0 HEAT=0 TANGENTIAL=0 INTRUDER=0
	$(MAKE) reconfigure-test THREE_D=1 STRESS=0 HEAT=0 TANGENTIAL=0 INTRUDER=0
	$(MAKE) reconfigure-test THREE_D=0 STRESS=1 HEAT=1 TANGENTIAL=0 INTRUDER=0
	$(MAKE) reconfigure-test THREE_D=1 STRESS=1 HEAT=1 TANGENTIAL=0 INTRUDER=0

graphics-performance-test: tests/test_graphics_performance.py
	$(PYTHON) -B $<

# Clean target
clean:
	rm -f a.out EDMD EDMD_debug EDMD_profile graphics_build *.o
	rm -rf $(GRAPHICS_BUILD_ROOT)
	rm -rf build/tests

-include $(GRAPHICS_DEPENDENCIES) $(WEB_DEPENDENCIES)


.DEFAULT_GOAL := EDMD
