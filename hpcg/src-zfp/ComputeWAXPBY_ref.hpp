
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

#ifndef COMPUTEWAXPBY_REF_HPP
#define COMPUTEWAXPBY_REF_HPP
#include "Vector.hpp"

#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This is the reference WAXPBY impmentation.  It CANNOT be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector.

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY
*/
template<class Datatype1, class Datatype2, class Datatype3>
int ComputeWAXPBY_ref(const local_int_t n, const double alpha,
    const Vector<Datatype1> & x, const double beta, const Vector<Datatype2> & y,
    Vector<Datatype3> & w) {

  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

  const zfp::array1<Datatype1>& xv = x.values;
  const zfp::array1<Datatype2>& yv = y.values;
  zfp::array1<Datatype3>& wv = w.values;

  if (alpha==1.0) {
//#ifndef HPCG_NO_OPENMP
//    #pragma omp parallel for
//#endif
    for (local_int_t i=0; i<n; i++) wv[i] = xv[i] + beta * yv[i];
  } else if (beta==1.0) {
//#ifndef HPCG_NO_OPENMP
//    #pragma omp parallel for
//#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + yv[i];
  } else  {
//#ifndef HPCG_NO_OPENMP
//    #pragma omp parallel for
//#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
  }

  return 0;
}

#endif // COMPUTEWAXPBY_REF_HPP
