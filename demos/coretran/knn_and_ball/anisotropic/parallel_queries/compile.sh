#!/bin/bash

PATH_KDTREE="../../../../../src/bin"
LIBS="-Wl,-rpath,$HOME/coretran/lib $HOME/coretran/lib/libcoretran.so"
INC="-I$HOME/coretran/library/include/coretran/"

gfortran -O3 -g -fopenmp -mieee-fp -ftree-vectorize -march=native -mcmodel=medium \
-J${PATH_KDTREE} ${INC} ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test ${LIBS}

