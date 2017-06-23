
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
 @file ComputeDotProduct.cpp

 HPCG routine
 */

#include "ComputeDotProduct.hpp"
#include "ComputeDotProduct_ref.hpp"

/*!
  Routine to compute the dot product of two vectors.

  This routine calls the reference dot-product implementation by default, but
  can be replaced by a custom routine that is optimized and better suited for
  the target system.

  @param[in]  n the number of vector elements (on this processor)
  @param[in]  x, y the input vectors
  @param[out] result a pointer to scalar value, on exit will contain the result.
  @param[out] time_allreduce the time it took to perform the communication between processes
  @param[out] isOptimized should be set to false if this routine uses the reference implementation (is not optimized); otherwise leave it unchanged

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct_ref
*/
template<class datatype1, class datatype2>
int ComputeDotProduct(const local_int_t n, const Vector<datatype1> & x,
    const Vector<datatype2> & y, double & result, double & time_allreduce,
    bool & isOptimized) {

  // This line and the next two lines should be removed and your version of ComputeDotProduct should be used.
  isOptimized = false;
  return ComputeDotProduct_ref(n, x, y, result, time_allreduce);
}



template int ComputeDotProduct<float, float>(const local_int_t n,
    const Vector<float> & x, const Vector<float> & y, double & result,
    double & time_allreduce, bool & isOptimized);
template int ComputeDotProduct<double, float>(const local_int_t n,
    const Vector<double> & x, const Vector<float> & y, double & result,
    double & time_allreduce, bool & isOptimized);
template int ComputeDotProduct<double, double>(const local_int_t n,
    const Vector<double> & x, const Vector<double> & y, double & result,
    double & time_allreduce, bool & isOptimized);
