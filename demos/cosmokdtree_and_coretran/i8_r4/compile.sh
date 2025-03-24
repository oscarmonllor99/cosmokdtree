#!/bin/bash
PATH_KDTREE="../../../src/bin"
LIBS_CORETRAN="/usr/local/lib/libcoretran.a -lcoretran"
INC_CORETRAN="-I/usr/local/include/coretran/"

#gfortran -O3 -g -fopenmp -J${PATH_KDTREE}  ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test
gfortran -O3 -g -fopenmp -J${PATH_KDTREE} -I${PATH_KDTREE} ${INC_CORETRAN} ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test ${LIBS_CORETRAN}