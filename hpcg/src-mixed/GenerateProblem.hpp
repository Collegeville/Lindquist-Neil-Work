
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

#ifndef GENERATEPROBLEM_HPP
#define GENERATEPROBLEM_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void GenerateProblem(SparseMatrix & A, Vector<float> * b, Vector<float> * x, Vector<float> * xexact);
#endif // GENERATEPROBLEM_HPP
