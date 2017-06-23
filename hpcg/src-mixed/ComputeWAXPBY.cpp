
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
 @file ComputeWAXPBY.cpp

 HPCG routine
 */

#include "ComputeWAXPBY.hpp"
#include "ComputeWAXPBY_ref.hpp"

/*!
  Routine to compute the update of a vector with the sum of two
  scaled vectors where: w = alpha*x + beta*y

  This routine calls the reference WAXPBY implementation by default, but
  can be replaced by a custom, optimized routine suited for
  the target system.

  @param[in] n the number of vector elements (on this processor)
  @param[in] alpha, beta the scalars applied to x and y respectively.
  @param[in] x, y the input vectors
  @param[out] w the output vector
  @param[out] isOptimized should be set to false if this routine uses the reference implementation (is not optimized); otherwise leave it unchanged

  @return returns 0 upon success and non-zero otherwise

  @see ComputeWAXPBY_ref
*/
template<class datatype1, class datatype2, class datatype3>
int ComputeWAXPBY(const local_int_t n, const double alpha,
    const Vector<datatype1> & x, const double beta, const Vector<datatype2> & y,
    Vector<datatype3> & w, bool & isOptimized) {

  // This line and the next two lines should be removed and your version of ComputeWAXPBY should be used.
  isOptimized = false;
  return ComputeWAXPBY_ref(n, alpha, x, beta, y, w);
}

template int ComputeWAXPBY<float, float, float>(const local_int_t n,
    const double alpha, const Vector<float> & x, const double beta,
    const Vector<float> & y, Vector<float> & w, bool & isOptimized);
    
template int ComputeWAXPBY<float, float, double>(const local_int_t n,
    const double alpha, const Vector<float> & x, const double beta,
    const Vector<float> & y, Vector<double> & w, bool & isOptimized);
    
template int ComputeWAXPBY<double, float, float>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<float> & y, Vector<float> & w, bool & isOptimized);
    
template int ComputeWAXPBY<double, float, double>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<float> & y, Vector<double> & w, bool & isOptimized);
    
template int ComputeWAXPBY<double, double, float>(const local_int_t n,
    const double alpha, const Vector<double> & x, const double beta,
    const Vector<double> & y, Vector<float> & w, bool & isOptimized);
