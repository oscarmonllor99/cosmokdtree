#!/bin/bash
PATH_KDTREE="Path/To/cosmokdtree.mod" ## should be inside src/bin

gfortran -O3 -g -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -c FoF.f90 -o fof.o 
gfortran -O3 -g -mcmodel=medium -fopenmp -mieee-fp -ftree-vectorize -march=native -J${PATH_KDTREE}  ${PATH_KDTREE}/cosmokdtree.o fof.o test.f90 -o test
