
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

#ifndef COMPUTEMG_HPP
#define COMPUTEMG_HPP
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "ComputeMG_ref.hpp"

/*!
  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On exit contains the result of the multigrid V-cycle with r as the RHS, x is the approximation to Ax = r.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeMG_ref
*/
template<class Datatype1, class Datatype2, class Datatype3>
int ComputeMG(const SparseMatrix<Datatype1>  & A, const Vector<Datatype2> & r, Vector<Datatype3> & x) {

  // This line and the next two lines should be removed and your version of ComputeSYMGS should be used.
  A.isMgOptimized = false;
  return ComputeMG_ref(A, r, x);
}

#endif // COMPUTEMG_HPP
