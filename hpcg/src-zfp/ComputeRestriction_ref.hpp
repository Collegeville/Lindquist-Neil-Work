
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

#ifndef COMPUTERESTRICTION_REF_HPP
#define COMPUTERESTRICTION_REF_HPP

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "Vector.hpp"
#include "SparseMatrix.hpp"

/*!
  Routine to compute the coarse residual vector.

  @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf, the fine grid matrix-vector product and mgData->rc the coarse residual vector.
  @param[in]    rf - Fine grid RHS.


  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
template<class DatatypeA, class DatatypeR>
int ComputeRestriction_ref(const SparseMatrix<DatatypeA> & A, const Vector<DatatypeR> & rf) {

  zfp::array1<DatatypeA>& Axfv = A.mgData->Axf->values;
  const zfp::array1<DatatypeR>& rfv = rf.values;
  zfp::array1<DatatypeA>& rcv = A.mgData->rc->values;
  local_int_t * f2c = A.mgData->f2cOperator;
  local_int_t nc = A.mgData->rc->localLength;

//#ifndef HPCG_NO_OPENMP
//#pragma omp parallel for
//#endif
  for (local_int_t i=0; i<nc; ++i) rcv[i] = rfv[f2c[i]] - Axfv[f2c[i]];

  return 0;
}

#endif // COMPUTERESTRICTION_REF_HPP