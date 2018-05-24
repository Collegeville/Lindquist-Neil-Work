
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

#ifndef COMPUTESPMV_HPP
#define COMPUTESPMV_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"
#include "ComputeSPMV_ref.hpp"

/*!
  Routine to compute sparse matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This routine calls the reference SpMV implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV_ref
*/
template<class DatatypeA, class DatatypeX, class DatatypeY>
int ComputeSPMV( const SparseMatrix<DatatypeA> & A, Vector<DatatypeX> & x, Vector<DatatypeY> & y) {

  // This line and the next two lines should be removed and your version of ComputeSPMV should be used.
  A.isSpmvOptimized = false;
  return ComputeSPMV_ref(A, x, y);
}

#endif  // COMPUTESPMV_HPP
