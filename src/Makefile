#directories
BINDIR := bin
NAME := cosmokdtree

#compiler flags
PERIODIC := $(if $(PERIODIC),$(PERIODIC),0)
$(info PERIODIC = $(PERIODIC))

LONGINT := $(if $(LONGINT),$(LONGINT),0)
$(info LONGINT = $(LONGINT))

DOUBLEPRECISION := $(if $(DOUBLEPRECISION),$(DOUBLEPRECISION),0)
$(info DOUBLEPRECISION = $(DOUBLEPRECISION))

DIMEN := $(if $(DIMEN),$(DIMEN),3)
$(info DIMEN = $(DIMEN))

PYTHON := $(if $(PYTHON),$(PYTHON),0)
$(info PYTHON = $(PYTHON))

#setting DOUBLEPRECISION=1 and LONGINT=1 if PYTHON=1
ifeq ($(PYTHON),1)
	DOUBLEPRECISION := 1
	LONGINT := 1
endif

# Fortran flags
FC = gfortran
FFLAGS = -cpp -O0 -g -mcmodel=medium -fopenmp -march=native -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
#FFLAGS = -cpp -O3 -g -fopenmp -mieee-fp -ftree-vectorize -funroll-loops -march=native -mcmodel=medium
FFLAGS += -Dperiodic=$(PERIODIC) -Dlongint=$(LONGINT) -Ddoubleprecision=$(DOUBLEPRECISION) -Ddimen=$(DIMEN) 

#sources and .o .mod
MOD_SRC = $(NAME).f90
MOD_OBJ = $(BINDIR)/$(NAME).o
MOD_MOD = $(BINDIR)/$(NAME).mod  

#python wrapper
PYMOD_NAME := pycosmokdtree
PYMOD_SRC  := cosmokdtree2py.f90

# 'all' so Make knows python should be built if python=1
.PHONY: all clean python

ifeq ($(PYTHON),1)
all: $(MOD_OBJ) python
else
all: $(MOD_OBJ)
endif

######TARGETS

# Compile the MAIN module
$(MOD_OBJ): $(MOD_SRC) | $(BINDIR)
	$(info Compiling $(MOD_SRC) to $(MOD_OBJ) ...)
	$(FC) $(FFLAGS) -J $(BINDIR) -c -o $(MOD_OBJ) $(MOD_SRC)

# Compile the Python wrapper using f2py
python: $(MOD_OBJ) $(PYMOD_SRC)
	$(info Building Python module $(PYMOD_NAME) ...)
	f2py --f90flags='$(FFLAGS) -I$(BINDIR)' -c -m $(PYMOD_NAME) $(MOD_OBJ) $(PYMOD_SRC)

# Ensure bin/ exists
$(BINDIR):
	mkdir -p $(BINDIR)

# Clean
clean:
	rm -f $(MOD_OBJ) $(MOD_MOD)
	rm -f $(PYMOD_NAME)*.so
	rm -rf build

info:
	@echo ""
	@echo "***********************************************************"
	@echo "***                    cosmokdtree                     ***"
	@echo "***********************************************************"
	@echo "* Óscar Monllor-Berbegal / Univ. de València, 2025  *"
	@echo "***********************************************************"
	@echo "*** FLAGS for compiling ***"
	@echo "- PERIODIC (optional, default: 0): set to 1 if you want to use periodic boundary conditions"
	@echo "- LONGINT (optional, default: 0): set to 1 if you want to use long integers (8 bytes) instead of integers (4 bytes)"
	@echo "- DOUBLEPRECISION (optional, default: 0): set to 1 if you want to use double precision (8 bytes)"
	@echo "- DIMEN (optional, default: 3): set to the number of dimensions of the data (2,3,4,...)"
	@echo "- PYTHON (optional, default: 0): set to 1 if you want to build the Python wrapper"
