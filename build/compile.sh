#!/bin/bash


gfortran -c ../External/f90getopt.F90 
gfortran -c ../src/header.f90 
gfortran -c ../src/parser.f90 
gfortran -c ../src/force_field_calc.f90 
gfortran -c ../src/main.f90 
gfortran *.o -o ../md_program &&
echo "All compiled"


