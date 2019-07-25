#!/bin/bash

gcc pdgesv_example.c -lscalapack-openmpi -lblas -lpthread -lgfortran
mpirun -np 6 --oversubscribe ./a.out 9 2 2 3
