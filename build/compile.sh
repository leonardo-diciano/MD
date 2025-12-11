#!/bin/bash

rm definitions.o definitions.mod f90getopt.o force_field_mod.mod  header.o minimization_mod.mod parser_mod.mod f90getopt.mod force_field_calc.o header_mod.mod main.o minimization.o parser.o

gfortran -c ../src/definitions.f90 
gfortran -c ../src/print_mod.f90 
gfortran -c ../src/lin_alg.f90 
gfortran -c ../External/f90getopt.F90 
gfortran -c ../src/header.f90 
gfortran -c ../src/parser.f90 
gfortran -c ../src/force_field_calc.f90 
gfortran -c ../src/minimization.f90 
gfortran -c ../src/main.f90 
gfortran *.o -o ../md_program &&
echo "All compiled"


