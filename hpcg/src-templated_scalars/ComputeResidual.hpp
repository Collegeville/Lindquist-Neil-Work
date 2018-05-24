
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

#ifndef COMPUTERESIDUAL_HPP
#define COMPUTERESIDUAL_HPP
#ifndef HPCG_NO_MPI
#include <mpi.h>
#endif
#ifndef HPCG_NO_OPENMP
#include <omp.h>
#endif

#include "Vector.hpp"

#ifdef HPCG_DETAILED_DEBUG
#include <fstream>
#include "hpcg.hpp"
#endif

#include <cmath>  // needed for fabs
#include "Vector.hpp"
#ifdef HPCG_DETAILED_DEBUG
#include <iostream>
#endif

/*!
  Routine to compute the inf-norm difference between two vectors where:

  @param[in]  n        number of vector elements (local to this processor)
  @param[in]  v1, v2   input vectors
  @param[out] residual pointer to scalar value; on exit, will contain result: inf-norm difference

  @return Returns zero on success and a non-zero value otherwise.
*/
template<class Datatype1, class Datatype2, class ResultDatatype>
int ComputeResidual(const local_int_t n, const Vector<Datatype1> & v1,
                    const Vector<Datatype2> & v2, ResultDatatype & residual) {

  Datatype1 * v1v = v1.values;
  Datatype2 * v2v = v2.values;
  ResultDatatype local_residual = 0.0;

#ifndef HPCG_NO_OPENMP
  #pragma omp parallel default(none) shared(local_residual, v1v, v2v)
  {
    ResultDatatype threadlocal_residual = 0.0;
    #pragma omp for
    for (local_int_t i=0; i<n; i++) {
      ResultDatatype diff = std::fabs(ResultDatatype(v1v[i]) - ResultDatatype(v2v[i]));
      if (diff > threadlocal_residual) threadlocal_residual = diff;
    }
    #pragma omp critical
    {
      if (threadlocal_residual>local_residual) local_residual = threadlocal_residual;
    }
  }
#else // No threading
  for (local_int_t i=0; i<n; i++) {
    ResultDatatype diff = std::fabs(ResultDatatype(v1v[i]) - ResultDatatype(v2v[i]));
    if (diff > local_residual) local_residual = diff;
#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout << " Computed, exact, diff = " << v1v[i] << " " << v2v[i] << " " << diff << std::endl;
#endif
  }
#endif

#ifndef HPCG_NO_MPI
  // Use MPI's reduce function to collect all partial sums
  const MPITypeWrapper mpitype = getMPIType(local_residual);
  ResultDatatype global_residual = 0;
  MPI_Allreduce(&local_residual, &global_residual, mpitype.count, mpitype.mpitype, MPI_MAX, MPI_COMM_WORLD);
  residual = global_residual;
#else
  residual = local_residual;
#endif

  return 0;
}

#endif // COMPUTERESIDUAL_HPP