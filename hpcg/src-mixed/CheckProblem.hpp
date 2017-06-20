
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

#ifndef CHECKPROBLEM_HPP
#define CHECKPROBLEM_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"

void CheckProblem(SparseMatrix<float> & A, Vector<float> * b, Vector<float> * x, Vector<float> * xexact);
#endif // CHECKPROBLEM_HPP
