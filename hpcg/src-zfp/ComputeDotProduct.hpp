
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

#ifndef COMPUTEDOTPRODUCT_HPP
#define COMPUTEDOTPRODUCT_HPP
#include "Vector.hpp"
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
template<class Datatype1, class Datatype2, class Datatype3>
int ComputeDotProduct(const local_int_t n, Vector<Datatype1> & x,
    Vector<Datatype2> & y, Datatype3 & result, double & time_allreduce,
    bool & isOptimized) {

  // This line and the next two lines should be removed and your version of ComputeDotProduct should be used.
  isOptimized = false;
  return ComputeDotProduct_ref(n, x, y, result, time_allreduce);
}

#endif // COMPUTEDOTPRODUCT_HPP
