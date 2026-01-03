#!/bin/bash

rm *.o *.mod

gfortran -c ../src/definitions.f90 
gfortran -c ../src/print_mod.f90 
gfortran -c ../src/lin_alg.f90 
gfortran -c ../External/f90getopt.F90 
gfortran -c ../src/force_field_calc.f90 
gfortran -c ../src/header.f90 
gfortran -c ../src/parser.f90 
gfortran -c ../src/pbc.f90 
gfortran -c ../src/minimization.f90 
gfortran -c ../src/propagators.f90 
gfortran -c ../src/ensemble.f90 
gfortran -c ../src/simulation_subroutines.f90 
gfortran -c ../src/metadynamics.f90 
gfortran -c ../src/simulation_mod.f90 
gfortran -c ../src/main.f90 
gfortran *.o -o ../md_program &&
echo "All compiled"

