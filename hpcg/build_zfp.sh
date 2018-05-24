#! /bin/sh

#terminate script on failure
set -e

echo
echo Starting Serial
echo
cd Linux_Serial_ZFP
touch src/garbage
rm src/*
make
cd ..


echo
echo Starting Intel MPI
echo
cd MPI_Intel_OMP_ZFP
touch src/garbage
rm src/*
make
cd ..

