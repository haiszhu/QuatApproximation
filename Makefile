# QuatApproximation Makefile
# Builds: libQuatApproximation.a + QuatApproximation_mex.<EXT>
#
# Parallel-independent sibling of LineQuaaadrature-legacy (no cross-package
# Fortran `use`, no cross-package linking; see the design spec at
# docs/superpowers/specs/2026-05-15-QuatApproximation-legacy-layout-design.md).
#
# Usage:
#   make          — show targets
#   make lib      — build build/libQuatApproximation.a
#   make mex      — build matlab/QuatApproximation_mex.<EXT>
#   make clean    — remove build artifacts

SHELL := /bin/bash

ROOT       := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ROOT       := $(patsubst %/,%,$(ROOT))
SRC_DIR    := $(ROOT)/src
MATLAB_DIR := $(ROOT)/matlab
BLD_DIR    := $(ROOT)/build

# ---- compilers ----
# Homebrew retires versioned formulae -- gcc@15 is now a dangling symlink into
# gcc 16 -- so a hardcoded gfortran-NN breaks on the next upgrade.  Take the
# newest one installed; FC/CC given on the command line still win.
ifneq ($(filter default environment,$(origin FC)),)
  FC := $(shell for c in gfortran-16 gfortran-15 gfortran-14 gfortran-13 gfortran; do \
                  command -v $$c >/dev/null 2>&1 && { echo $$c; exit 0; }; done; echo gfortran)
endif
ifneq ($(filter default environment,$(origin CC)),)
  CC := $(shell for c in gcc-16 gcc-15 gcc-14 gcc-13 gcc; do \
                  command -v $$c >/dev/null 2>&1 && { echo $$c; exit 0; }; done; echo gcc)
endif
MW := ~/mwrap/mwrap
MWFLAGS := -c99complex -i8 -mex

# ---- platform detection ----
UNAME := $(shell uname)
ARCH  := $(shell uname -m)

# Homebrew prefix: /opt/homebrew on Apple Silicon, /usr/local on Intel.
BREW := /opt/homebrew
ifeq ($(UNAME), Darwin)
  ifneq ($(ARCH), arm64)
    BREW := /usr/local
  endif
endif

ifeq ($(UNAME), Darwin)
  MATLAB_ROOT  := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  ifeq ($(ARCH), arm64)
    MATLAB_ARCH  := maca64
    MEX_EXT      := mexmaca64
  else
    MATLAB_ARCH  := maci64
    MEX_EXT      := mexmaci64
  endif
  OPENBLAS_DIR := $(BREW)/opt/openblas-singlethread
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  MEX_LDFLAGS  := -bundle -Wl,-undefined,dynamic_lookup
else
  MATLAB_ROOT  := $(shell ls -d /usr/local/MATLAB/R* 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH  := glnxa64
  MEX_EXT      := mexa64
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmx -lmex -lmat -lm
  MEX_LDFLAGS  := -shared
endif

OPENBLAS_LIBS := -L$(OPENBLAS_DIR)/lib -lopenblas
HDF5_ROOT ?= $(BREW)/opt/hdf5
HDF5_INC := -I$(HDF5_ROOT)/include
HDF5_LIBS := -L$(HDF5_ROOT)/lib -lhdf5

# ---- Fortran flags (identical to LineQuaaadrature) ----
# -fdefault-integer-8 : 8-byte integers (matches mwrap -i8)
FFLAGS := -g -O3 -fPIC \
           -ffp-contract=off -fno-unsafe-math-optimizations \
           -fallow-argument-mismatch -std=legacy -w \
           -fdefault-integer-8 \
           -frecursive \
           -cpp \
           -march=native -funroll-loops \
           -fopenmp \
           -J$(BLD_DIR) -I$(BLD_DIR) -I$(SRC_DIR)

# ---- sources and objects ----
# Add new _mod.f90 / _mex.f90 pairs here as routines are ported.
QA_SOURCES := $(SRC_DIR)/quatapproximation_mod.f90 \
              $(SRC_DIR)/quatapproximation_mex.f90 \
              $(SRC_DIR)/harmonic_mod.f90 \
              $(SRC_DIR)/harmonic_mex.f90 \
              $(SRC_DIR)/koorn_geom_mod.f90 \
              $(SRC_DIR)/koorn_geom_mex.f90 \
              $(SRC_DIR)/qkernel_mod.f90 \
              $(SRC_DIR)/qkernel_mex.f90 \
              $(SRC_DIR)/omega_mod.f90 \
              $(SRC_DIR)/omega_mex.f90 \
              $(SRC_DIR)/tensor_geom_mod.f90 \
              $(SRC_DIR)/tensor_geom_mex.f90

QA_OBJECTS := $(patsubst $(SRC_DIR)/%.f90, $(BLD_DIR)/%.o, $(QA_SOURCES)) \
              $(BLD_DIR)/hdf5_io.o

LIB := $(BLD_DIR)/libQuatApproximation.a

# ---- mwrap-generated gateway ----
MW_SRC  := $(MATLAB_DIR)/QuatApproximation.mw
MEX_C   := $(MATLAB_DIR)/QuatApproximation_mex.c
MEX_OUT := $(MATLAB_DIR)/QuatApproximation_mex.$(MEX_EXT)

# ---- test binaries ----
TEST_PARABOLOID_COND_SRC := $(ROOT)/test/harmonic_approx/test_paraboloid_cond.f90
TEST_PARABOLOID_COND_BIN := $(BLD_DIR)/test_paraboloid_cond
TEST_PARABOLOID_REFINE_SRC := $(ROOT)/test/harmonic_approx/test_paraboloid_refinement.f90
TEST_PARABOLOID_REFINE_BIN := $(BLD_DIR)/test_paraboloid_refinement

# ============================================================

.PHONY: all mex lib test_paraboloid_cond test_paraboloid_refinement clean

all:
	@echo "QuatApproximation local build targets"
	@echo ""
	@echo "  make mex                          build $(notdir $(MEX_OUT))"
	@echo "  make lib                          build $(notdir $(LIB))"
	@echo "  make test_paraboloid_cond         build build/test_paraboloid_cond"
	@echo "  make test_paraboloid_refinement   build build/test_paraboloid_refinement"
	@echo "  make clean                        remove build artifacts"
	@echo ""

mex: $(MEX_OUT)

lib: $(LIB)

test_paraboloid_cond: $(TEST_PARABOLOID_COND_BIN)

$(TEST_PARABOLOID_COND_BIN): $(LIB) $(TEST_PARABOLOID_COND_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST_PARABOLOID_COND_SRC) \
	  -L$(BLD_DIR) -lQuatApproximation \
	  -framework Accelerate -lgfortran -lquadmath -lm -o $(TEST_PARABOLOID_COND_BIN)

test_paraboloid_refinement: $(TEST_PARABOLOID_REFINE_BIN)

$(TEST_PARABOLOID_REFINE_BIN): $(LIB) $(TEST_PARABOLOID_REFINE_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST_PARABOLOID_REFINE_SRC) \
	  -L$(BLD_DIR) -lQuatApproximation \
	  -framework Accelerate -lgfortran -lquadmath -lm -o $(TEST_PARABOLOID_REFINE_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# ---- compile Fortran objects ----
$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BLD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# ---- compile hdf5_io.c (design spec §13) ----
# Standalone build (default): exported names are bare hdf5_* — Fortran
# bind(C) interfaces use plain hdf5_<name> with no name= clause.
# External-composition build (EXTERNAL_LINK=1): exported names become
# qa_hdf5_* via -DHDF5_PREFIX=qa_, so libQuatApproximation.a can be
# linked alongside libLineQuaaadrature.a without symbol collision. In
# that mode, Fortran bind(C) interfaces need an explicit
# name='qa_hdf5_<name>' (one-time source patch, deferred until needed).
ifeq ($(EXTERNAL_LINK),1)
  HDF5_PREFIX_FLAG := -DHDF5_PREFIX=qa_
else
  HDF5_PREFIX_FLAG :=
endif

$(BLD_DIR)/hdf5_io.o: $(SRC_DIR)/hdf5_io.c | $(BLD_DIR)
	$(CC) -c -fPIC $(HDF5_PREFIX_FLAG) $(HDF5_INC) $< -o $@

# ---- module dependency graph ----
$(BLD_DIR)/quatapproximation_mex.o:  $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/harmonic_mex.o:           $(BLD_DIR)/harmonic_mod.o
$(BLD_DIR)/koorn_geom_mex.o:         $(BLD_DIR)/koorn_geom_mod.o $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/qkernel_mod.o:            $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/qkernel_mex.o:            $(BLD_DIR)/qkernel_mod.o $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/omega_mod.o:              $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/omega_mex.o:              $(BLD_DIR)/omega_mod.o $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/tensor_geom_mod.o:        $(BLD_DIR)/quatapproximation_mod.o
$(BLD_DIR)/tensor_geom_mex.o:        $(BLD_DIR)/tensor_geom_mod.o $(BLD_DIR)/quatapproximation_mod.o

# ---- static library ----
$(LIB): $(QA_OBJECTS)
	ar rcs $@ $^

# ---- mwrap: two-pass generation ----
$(MEX_C): $(MW_SRC) | $(BLD_DIR)
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) QuatApproximation_mex -mb -list QuatApproximation.mw
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) QuatApproximation_mex -c QuatApproximation_mex.c QuatApproximation.mw
	perl -pi -e 's/_{2,}/_/g' $(MEX_C)

# ---- link MEX (via gcc-15 directly, not MATLAB's `mex` script) ----
$(MEX_OUT): $(LIB) $(MEX_C)
	$(CC) $(MEX_LDFLAGS) -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  $(MATLAB_INC) \
	  $(MEX_C) \
	  -L$(BLD_DIR) -lQuatApproximation \
	  $(HDF5_LIBS) \
	  $(MATLAB_LIBS) \
	  $(OPENBLAS_LIBS) \
	  -lgfortran -lquadmath -lm \
	  -o $(MEX_OUT)

# ---- clean ----
clean:
	rm -rf $(BLD_DIR) $(MEX_OUT)
