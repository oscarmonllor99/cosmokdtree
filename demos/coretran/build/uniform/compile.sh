#!/bin/bash

LIBS="-Wl,-rpath,$HOME/coretran/lib $HOME/coretran/lib/libcoretran.so"
INC="-I$HOME/coretran/library/include/coretran/"

gfortran -O3 -g -mieee-fp -ftree-vectorize -march=native -mcmodel=medium ${INC} tree_test.f90 -o tree_test ${LIBS}

