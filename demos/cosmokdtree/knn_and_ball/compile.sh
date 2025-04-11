#!/bin/bash
PATH_KDTREE="../../../src/bin"

gfortran -O3 -g -fopenmp -J${PATH_KDTREE}  ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test
## gfortran -O1 -g -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -fcheck=all -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -J${PATH_KDTREE}  ${PATH_KDTREE}/cosmokdtree.o tree_test.f90 -o tree_test
