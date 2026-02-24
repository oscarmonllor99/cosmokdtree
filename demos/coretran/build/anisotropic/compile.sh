#!/bin/bash

LIBS="-Wl,-rpath,$HOME/coretran/lib $HOME/coretran/lib/libcoretran.so"
INC="-I$HOME/coretran/library/include/coretran/"

gfortran -O3 -g -fopenmp -mieee-fp -ftree-vectorize -march=native -mcmodel=medium \
 ${INC} tree_test.f90 -o tree_test ${LIBS}

##gfortran -O3 -g -fopenmp -mieee-fp -ftree-vectorize -march=native -mcmodel=medium \
## -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -J${PATH_KDTREE} ${INC} ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test ${LIBS}
