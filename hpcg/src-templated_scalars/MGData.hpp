
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
 @file MGData.hpp

 HPCG data structure
 */

#ifndef MGDATA_HPP
#define MGDATA_HPP

#include <cassert>
#include "SparseMatrix.hpp"
#include "Vector.hpp"

template<class ResidualDatatype, class SolutionDatatype>
struct MGData {
  int numberOfPresmootherSteps; // Call ComputeSYMGS this many times prior to coarsening
  int numberOfPostsmootherSteps; // Call ComputeSYMGS this many times after coarsening
  local_int_t * f2cOperator; //!< 1D array containing the fine operator local IDs that will be injected into coarse space.
  Vector<ResidualDatatype> * rc; // coarse grid residual vector
  Vector<SolutionDatatype> * xc; // coarse grid solution vector
  Vector<ResidualDatatype> * Axf; // fine grid residual vector
  /*!
   This is for storing optimized data structres created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  void * optimizationData;
};

/*!
 Constructor for the data structure of CG vectors.

 @param[in] Ac - Fully-formed coarse matrix
 @param[in] f2cOperator -
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
template<class ResidualDatatype, class SolutionDatatype>
inline void InitializeMGData(local_int_t * f2cOperator, Vector<ResidualDatatype> * rc, Vector<SolutionDatatype> * xc, Vector<ResidualDatatype> * Axf, MGData<ResidualDatatype, SolutionDatatype> & data) {
  data.numberOfPresmootherSteps = 1;
  data.numberOfPostsmootherSteps = 1;
  data.f2cOperator = f2cOperator; // Space for injection operator
  data.rc = rc;
  data.xc = xc;
  data.Axf = Axf;
  return;
}

/*!
 Destructor for the CG vectors data.

 @param[inout] data the MG data structure whose storage is deallocated
 */
template<class ResidualDatatype, class SolutionDatatype>
inline void DeleteMGData(MGData<ResidualDatatype, SolutionDatatype> & data) {

  delete [] data.f2cOperator;
  DeleteVector(*data.Axf);
  DeleteVector(*data.rc);
  DeleteVector(*data.xc);
  delete data.Axf;
  delete data.rc;
  delete data.xc;
  return;
}

#endif // MGDATA_HPP

