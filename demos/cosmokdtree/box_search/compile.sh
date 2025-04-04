#!/bin/bash
PATH_KDTREE="../../../src/bin"

gfortran -O3 -g -fopenmp -J${PATH_KDTREE}  ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test
