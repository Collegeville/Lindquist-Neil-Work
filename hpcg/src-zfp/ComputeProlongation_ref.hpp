
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

#ifndef COMPUTEPROLONGATION_REF_HPP
#define COMPUTEPROLONGATION_REF_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif


/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
template<class DatatypeA, class DatatypeX>
int ComputeProlongation_ref(const SparseMatrix<DatatypeA> & Af, Vector<DatatypeX> & xf) {

  zfp::array1<DatatypeX>& xfv = xf.values;
  zfp::array1<DatatypeA>& xcv = Af.mgData->xc->values;
  local_int_t * f2c = Af.mgData->f2cOperator;
  local_int_t nc = Af.mgData->rc->localLength;

//#ifndef HPCG_NO_OPENMP
//#pragma omp parallel for
//#endif
// TODO: Somehow note that this loop can be safely vectorized since f2c has no repeated indices
  for (local_int_t i=0; i<nc; ++i) xfv[f2c[i]] += xcv[i]; // This loop is safe to vectorize

  return 0;
}

#endif // COMPUTEPROLONGATION_REF_HPP
