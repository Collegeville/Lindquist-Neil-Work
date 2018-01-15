#! /bin/tcsh

#clean old files
rm -fr CMakeFiles CMakeCache.txt cmake_install.cmake Makefile petra_power_method_LL

#run cmake
cmake -DTrilinos_DIR=~/cs_research/Trilinos/lib/cmake/Trilinos/ src/


#run make
make
