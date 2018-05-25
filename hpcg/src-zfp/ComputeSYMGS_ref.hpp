
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

#ifndef HPCG_NO_MPI
#include "ExchangeHalo.hpp"
#endif
#include <cassert>

/*!
  Computes one step of symmetric Gauss-Seidel:

  Assumption about the structure of matrix A:
  - Each row 'i' of the matrix has nonzero diagonal value whose address is matrixDiagonal[i]
  - Entries in row 'i' are ordered such that:
       - lower triangular terms are stored before the diagonal element.
       - upper triangular terms are stored after the diagonal element.
       - No other assumptions are made about entry ordering.

  Symmetric Gauss-Seidel notes:
  - We use the input vector x as the RHS and start with an initial guess for y of all zeros.
  - We perform one forward sweep.  x should be initially zero on the first GS sweep, but we do not attempt to exploit this fact.
  - We then perform one back sweep.
  - For simplicity we include the diagonal contribution in the for-j loop, then correct the sum after

  @param[in] A the known system matrix
  @param[in] r the input vector
  @param[inout] x On entry, x should contain relevant values, on exit x contains the result of one symmetric GS sweep with r as the RHS.


  @warning Early versions of this kernel (Version 1.1 and earlier) had the r and x arguments in reverse order, and out of sync with other kernels.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeSYMGS
*/
template<class DatatypeA, class DatatypeR, class DatatypeX>
int ComputeSYMGS_ref(const SparseMatrix<DatatypeA> & A, const Vector<DatatypeR> & r, Vector<DatatypeX> & x) {

  assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

#ifndef HPCG_NO_MPI
  ExchangeHalo(A,x);
#endif

  const local_int_t nrow = A.localNumberOfRows;
  local_int_t nx = A.geom->nx;
  local_int_t ny = A.geom->ny;
  DatatypeA ** matrixDiagonal = A.matrixDiagonal;  // An array of pointers to the diagonal entries A.matrixValues
  const zfp::array3<DatatypeR>& rv = r.values;
  zfp::array3<DatatypeX>& xv = x.values;

  for (local_int_t j=0; j< nrow; j++)  {
    local_int_t block = j>>6;
    local_int_t iz = (j>>4)&(local_int_t)0x03 + (block>>2)&(local_int_t(-1-15));
    local_int_t iy = (j>>2)&(local_int_t)0x03 + block&(local_int_t(-1-15));
    local_int_t ix = j&(local_int_t)0x03 + block<<2;
    local_int_t i = ix + iy*nx+iz*nx+ny;

    const DatatypeA * const currentValues = A.matrixValues[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const DatatypeA  currentDiagonal = matrixDiagonal[i][0]; // Current diagonal value
    DatatypeX sum = rv[i]; // RHS value

    for (int j=0; j< currentNumberOfNonzeros; j++) {
      local_int_t curCol = currentColIndices[j];
      sum -= DatatypeX(currentValues[j]) * xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;

  }

  // Now the back sweep.

  for (local_int_t j=0; j< nrow; j++)  {
    local_int_t block = j>>6;
    local_int_t iz = (j>>4)&(local_int_t)0x03 + (block>>2)&(local_int_t(-1-15));
    local_int_t iy = (j>>2)&(local_int_t)0x03 + block&(local_int_t(-1-15));
    local_int_t ix = j&(local_int_t)0x03 + block<<2;
    local_int_t i = ix + iy*nx+iz*nx+ny;

    const DatatypeA * const currentValues = A.matrixValues[i];
    const local_int_t * const currentColIndices = A.mtxIndL[i];
    const int currentNumberOfNonzeros = A.nonzerosInRow[i];
    const DatatypeA  currentDiagonal = matrixDiagonal[i][0]; // Current diagonal value
    DatatypeX sum = rv[i]; // RHS value

    for (int j = 0; j< currentNumberOfNonzeros; j++) {
      local_int_t curCol = currentColIndices[j];
      sum -= DatatypeX(currentValues[j])*xv[curCol];
    }
    sum += xv[i]*currentDiagonal; // Remove diagonal contribution from previous loop

    xv[i] = sum/currentDiagonal;
  }

  return 0;
}

#endif // COMPUTESYMGS_REF_HPP
