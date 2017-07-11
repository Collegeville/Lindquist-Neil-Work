
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

#ifndef COMPUTESYMGS_REF_HPP
#define COMPUTESYMGS_REF_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

extern double inner_time;
extern double first_inner_time;
extern double second_inner_time;
extern int calls;

int ComputeSYMGS_ref( const SparseMatrix  & A, const Vector & r, Vector & x);

#endif // COMPUTESYMGS_REF_HPP
