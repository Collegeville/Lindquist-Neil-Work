
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

#ifndef COMPUTESPMV_REF_HPP
#define COMPUTESPMV_REF_HPP
#include "Vector.hpp"
#include "SparseMatrix.hpp"

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>

/*!
  Routine to compute matrix vector product y = Ax where:
  Precondition: First call exchange_externals to get off-processor values of x

  This is the reference SPMV implementation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in]  A the known system matrix
  @param[in]  x the known vector
  @param[out] y the On exit contains the result: Ax.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSPMV
*/
template<class DatatypeA, class DatatypeX, class DatatypeY>
int ComputeSPMV_ref( const SparseMatrix<DatatypeA> & A, Vector<DatatypeX> & x, Vector<DatatypeY> & y) {

  assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
  assert(y.localLength>=A.localNumberOfRows);

#ifndef HPCG_NO_MPI
    ExchangeHalo(A,x);
#endif
  const zfp::array3<DatatypeX>& xv = x.values;
  zfp::array3<DatatypeY>& yv = y.values;
  const local_int_t nrow = A.localNumberOfRows;
//#ifndef HPCG_NO_OPENMP
//  #pragma omp parallel for
//#endif
  for (local_int_t i=0; i< nrow; i++)  {
    DatatypeY sum = 0.0;
    const DatatypeA * const cur_vals = A.matrixValues[i];
    const local_int_t * const cur_inds = A.mtxIndL[i];
    const int cur_nnz = A.nonzerosInRow[i];

    for (int j=0; j< cur_nnz; j++)
      sum += DatatypeY(cur_vals[j])*DatatypeY(xv[cur_inds[j]]);
    yv[i] = sum;
  }
  return 0;
}

#endif  // COMPUTESPMV_REF_HPP
