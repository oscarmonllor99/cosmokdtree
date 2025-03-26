# Directories and files
BINDIR := bin
NAME := cosmokdtree

# Compiler and flags
FC = gfortran
FFLAGS = -O3 -g -fopenmp -J $(BINDIR)  # -J bin ensures .mod files are stored in bin/

# Module source file
MOD_SRC = $(NAME).f90
MOD_OBJ = $(BINDIR)/$(NAME).o
MOD_MOD = $(BINDIR)/$(NAME).mod  # This will be automatically created by the compiler

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
