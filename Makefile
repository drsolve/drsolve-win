.DEFAULT_GOAL := all

CC = x86_64-w64-mingw32-gcc
AR = x86_64-w64-mingw32-ar

SRC_DIR := src
INCLUDE_DIR := include
WIN_GUI_DIR := win_gui
BUILD_DIR ?= build
THIRD_PARTY_DIR := third_party
BIN_DIR ?= bin
DLL_DIR ?= dll
LIB_DIR ?= lib

MINGW_INCLUDE_DIR := $(THIRD_PARTY_DIR)/mingw/include
FLINT_LIB_DIR := $(THIRD_PARTY_DIR)/lib
PML_INCLUDE_DIR := $(THIRD_PARTY_DIR)/pml/include
PML_LIB_DIR := $(THIRD_PARTY_DIR)/lib

DIXON_TARGET ?= dixon.exe
DIXON_REAL_TARGET ?= $(BIN_DIR)/dixon_cli_real.exe
DIXON_STATIC_LIB ?= $(LIB_DIR)/libdixon-1.a
DIXON_DLL ?= $(DLL_DIR)/libdixon-1.dll
DIXON_DLL_IMPORT ?= $(LIB_DIR)/libdixon-1.dll.a
WIN_GUI_TARGET ?= dixon_win_gui.exe
LAUNCHER_SOURCE := $(WIN_GUI_DIR)/dixon_launcher.c

COMMON_CFLAGS ?= -O2 -fopenmp -msse4.1 -mpclmul
GUI_CFLAGS ?= -O2 -D_WIN32_WINNT=0x0601
CPPFLAGS += -DHAVE_FLINT -DHAVE_PML
INCLUDE_FLAGS := -I$(INCLUDE_DIR) -I$(PML_INCLUDE_DIR) -I$(MINGW_INCLUDE_DIR)
CLI_LIB_DIRS := -L$(FLINT_LIB_DIR) -L$(PML_LIB_DIR)
CLI_LIBS := -lflint -lpml -lmpfr -lgmp -lopenblas -lm -lpthread
GUI_LIBS := -lcomctl32 -lcomdlg32 -lshlwapi -lwinmm

MATH_SOURCES = $(SRC_DIR)/dixon_complexity.c \
               $(SRC_DIR)/dixon_flint.c \
               $(SRC_DIR)/dixon_interface_flint.c \
               $(SRC_DIR)/dixon_test.c \
               $(SRC_DIR)/dixon_with_ideal_reduction.c \
               $(SRC_DIR)/fq_mat_det.c \
               $(SRC_DIR)/fq_mpoly_mat_det.c \
               $(SRC_DIR)/fq_multivariate_interpolation.c \
               $(SRC_DIR)/fq_mvpoly.c \
               $(SRC_DIR)/fq_nmod_roots.c \
               $(SRC_DIR)/fq_poly_mat_det.c \
               $(SRC_DIR)/fq_sparse_interpolation.c \
               $(SRC_DIR)/fq_unified_interface.c \
               $(SRC_DIR)/gf2n_mpoly.c \
               $(SRC_DIR)/gf2n_field.c \
               $(SRC_DIR)/gf2n_poly.c \
               $(SRC_DIR)/polynomial_system_solver.c \
               $(SRC_DIR)/resultant_with_ideal_reduction.c \
               $(SRC_DIR)/unified_mpoly_det.c \
               $(SRC_DIR)/unified_mpoly_interface.c \
               $(SRC_DIR)/unified_mpoly_resultant.c \
               $(SRC_DIR)/rational_system_solver.c \
               $(SRC_DIR)/fmpq_acb_roots.c

MATH_OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(MATH_SOURCES))
WIN_GUI_SOURCES := $(WIN_GUI_DIR)/dixon_gui.c $(WIN_GUI_DIR)/dixon_wrapper.c

.PHONY: all clean show-config package-layout

all: $(DIXON_TARGET) $(DIXON_REAL_TARGET) $(WIN_GUI_TARGET) $(DIXON_DLL)

$(BUILD_DIR) $(BIN_DIR) $(DLL_DIR) $(LIB_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(COMMON_CFLAGS) $(CPPFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

$(DIXON_STATIC_LIB): $(MATH_OBJECTS) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(DIXON_DLL): $(MATH_OBJECTS) | $(DLL_DIR) $(LIB_DIR)
	$(CC) -shared -Wl,--out-implib,$(DIXON_DLL_IMPORT) -o $@ $^ $(CLI_LIB_DIRS) $(CLI_LIBS) -fopenmp

$(DIXON_REAL_TARGET): dixon.c $(DIXON_STATIC_LIB) | $(BIN_DIR)
	$(CC) $(COMMON_CFLAGS) $(CPPFLAGS) $(INCLUDE_FLAGS) -o $@ $< $(DIXON_STATIC_LIB) $(CLI_LIB_DIRS) $(CLI_LIBS) -fopenmp

$(DIXON_TARGET): $(LAUNCHER_SOURCE)
	$(CC) $(GUI_CFLAGS) -o $@ $<

$(WIN_GUI_TARGET): $(WIN_GUI_SOURCES)
	$(CC) $(GUI_CFLAGS) -mwindows -o $@ $(WIN_GUI_SOURCES) $(GUI_LIBS)

package-layout:
	@echo "root executables: $(DIXON_TARGET) $(WIN_GUI_TARGET)"
	@echo "internal cli     : $(DIXON_REAL_TARGET)"
	@echo "dlls             : $(DLL_DIR)/"
	@echo "link libraries   : $(LIB_DIR)/"

show-config:
	@echo "CC                = $(CC)"
	@echo "MINGW_INCLUDE_DIR = $(MINGW_INCLUDE_DIR)"
	@echo "FLINT_LIB_DIR     = $(FLINT_LIB_DIR)"
	@echo "PML_INCLUDE_DIR   = $(PML_INCLUDE_DIR)"
	@echo "PML_LIB_DIR       = $(PML_LIB_DIR)"
	@echo "DIXON_TARGET      = $(DIXON_TARGET)"
	@echo "DIXON_REAL_TARGET = $(DIXON_REAL_TARGET)"
	@echo "WIN_GUI_TARGET    = $(WIN_GUI_TARGET)"
	@echo "DLL_DIR           = $(DLL_DIR)"
	@echo "LIB_DIR           = $(LIB_DIR)"

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(LIB_DIR) $(DIXON_TARGET) $(WIN_GUI_TARGET) dixon_cli_real.exe
	rm -f libdixon-1.a libdixon-1.dll libdixon-1.dll.a libdixon-win.a
