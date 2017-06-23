
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

/*!
 @file ComputeWAXPBY_ref.cpp

 HPCG routine
 */

#include "ComputeWAXPBY_ref.hpp"
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
template<class datatype1, class datatype2, class datatype3>
int ComputeWAXPBY_ref(const local_int_t n, const double alpha,
    const Vector<datatype1> & x, const double beta, const Vector<datatype2> & y,
    Vector<datatype3> & w) {

  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

  const datatype1 * const xv = x.values;
  const datatype2 * const yv = y.values;
  datatype3 * const wv = w.values;

  if (alpha==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = xv[i] + beta * yv[i];
  } else if (beta==1.0) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + yv[i];
  } else  {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for
#endif
    for (local_int_t i=0; i<n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
  }

  return 0;
}


template int ComputeWAXPBY_ref<float, float, float>(const local_int_t n,
    const double alpha, const Vector<float> & x, const double beta,
    const Vector<float> & y, Vector<float> & w);
    
template int ComputeWAXPBY_ref<float, float, double>(const local_int_t n,
    const double alpha, const Vector<float> & x, const double beta,
    const Vector<float> & y, Vector<double> & w);
    
template int ComputeWAXPBY_ref<double, float, float>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<float> & y, Vector<float> & w);
    
template int ComputeWAXPBY_ref<double, float, double>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<float> & y, Vector<double> & w);
    
template int ComputeWAXPBY_ref<double, double, float>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<double> & y, Vector<float> & w);
