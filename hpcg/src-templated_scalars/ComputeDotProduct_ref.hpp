
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

#ifndef COMPUTEDOTPRODUCT_REF_HPP
#define COMPUTEDOTPRODUCT_REF_HPP
#ifndef HPCG_NO_MPI
#include <mpi.h>
#include "MPITypeLookup.hpp"
#include "mytimer.hpp"
#endif
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include "Vector.hpp"

/*!
  Routine to compute the dot product of two vectors where:

  This is the reference dot-product implementation.  It _CANNOT_ be modified for the
  purposes of this benchmark.

  @param[in] n the number of vector elements (on this processor)
  @param[in] x, y the input vectors
  @param[in] result a pointer to scalar value, on exit will contain result.
  @param[out] time_allreduce the time it took to perform the communication between processes

  @return returns 0 upon success and non-zero otherwise

  @see ComputeDotProduct
*/
template<class Datatype1, class Datatype2, class ResultDatatype>
int ComputeDotProduct_ref(const local_int_t n, const Vector<Datatype1> & x,
    const Vector<Datatype2> & y, ResultDatatype & result, double & time_allreduce) {
  assert(x.localLength>=n); // Test vector lengths
  assert(y.localLength>=n);

  ResultDatatype local_result = 0.0;
  Datatype1 * xv = x.values;
  Datatype2 * yv = y.values;
  if ((void*)yv==(void*)xv) {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for reduction (+:local_result)
#endif
    for (local_int_t i=0; i<n; i++) local_result += ResultDatatype(xv[i])*ResultDatatype(xv[i]);
  } else {
#ifndef HPCG_NO_OPENMP
    #pragma omp parallel for reduction (+:local_result)
#endif
    for (local_int_t i=0; i<n; i++) local_result += ResultDatatype(xv[i])*ResultDatatype(yv[i]);
  }

#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to collect all partial sums
  double t0 = mytimer();
  const MPITypeWrapper mpitype = getMPIType(local_result);
  ResultDatatype global_result = 0.0;
  MPI_Allreduce(&local_result, &global_result, mpitype.count, mpitype.mpitype, MPI_SUM,
      MPI_COMM_WORLD);
  result = global_result;
  time_allreduce += mytimer() - t0;
#else
  time_allreduce += 0.0;
  result = local_result;
#endif

  return 0;
}

#endif // COMPUTEDOTPRODUCT_REF_HPP
