#! /bin/bash

# This launches TestStarter.jl with MPI and the given number of processes

mpirun --mca pml ob1 --mca btl tcp,self -np $1 julia --color=yes -O3 TestStarter.jl "${@:2}"

