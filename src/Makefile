# Directories and files
BINDIR := bin
NAME := cosmokdtree

PERIODIC := $(if $(PERIODIC),$(PERIODIC),0)
$(info PERIODIC = $(PERIODIC))

LONGINT := $(if $(LONGINT),$(LONGINT),0)
$(info LONGINT = $(LONGINT))

DOUBLEPRECISION := $(if $(DOUBLEPRECISION),$(DOUBLEPRECISION),0)
$(info DOUBLEPRECISION = $(DOUBLEPRECISION))

DIMEN := $(if $(DIMEN),$(DIMEN),3)
$(info NDIM = $(DIMEN))

# Compiler and flags
FC = gfortran
FFLAGS = -cpp -O3 -g -fopenmp -mieee-fp -ftree-vectorize -march=native -mcmodel=medium -J $(BINDIR)  # -J bin ensures .mod files are stored in bin/

FFLAGS += -Dperiodic=$(PERIODIC) -Dlongint=$(LONGINT) -Ddoubleprecision=$(DOUBLEPRECISION) -Ddimen=$(DIMEN) 

# Module source file
MOD_SRC = $(NAME).f90
MOD_OBJ = $(BINDIR)/$(NAME).o
MOD_MOD = $(BINDIR)/$(NAME).mod  

# Compile the module (MOD_OBJ depends on MOD_SRC)
$(MOD_OBJ): $(MOD_SRC) $(BINDIR)
	$(info Compiling $(MOD_SRC) to $(MOD_OBJ) ...)
	$(FC) $(FFLAGS) -c -o $(MOD_OBJ) $(MOD_SRC)

# Ensure bin/ exists before compilation
$(BINDIR):
	mkdir -p $(BINDIR)

# Phony targets
.PHONY: all clean
all: $(MOD_OBJ)  # Ensures MOD_OBJ gets compiled

# Clean
clean:
	rm -f $(MOD_OBJ) $(MOD_MOD)

info:
	@echo ""
	@echo "***********************************************************"
	@echo "***                    cosmokdtree                     ***"
	@echo "***********************************************************"
	@echo "* Óscar Monllor-Berbegal et al., Univ. de València, 2025  *"
	@echo "***********************************************************"
	@echo "*** FLAGS for compiling ***"
	@echo "- PERIODIC (optional, default: 0): set to 1 if you want to use periodic boundary conditions"
	@echo "- LONGINT (optional, default: 0): set to 1 if you want to use long integers (8 bytes) instead of integers (4 bytes)"
	@echo "- DOUBLEPRECISION (optional, default: 0): set to 1 if you want to use double precision (8 bytes)"
	@echo "- DIMEN (optional, default: 3): set to the number of dimensions of the data (2,3,4,...)"
