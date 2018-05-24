
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
 @file TestCG.hpp

 HPCG data structure
 */

#ifndef TESTCG_HPP
#define TESTCG_HPP

#include "hpcg.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "CGData.hpp"
#include <fstream>
#include <iostream>
using std::endl;
#include <vector>
#include "hpcg.hpp"

#include "CG.hpp"


template<class Datatype>
struct TestCGData {
  int count_pass; //!< number of succesful tests
  int count_fail;  //!< number of succesful tests
  int expected_niters_no_prec; //!< expected number of test CG iterations without preconditioning with diagonally dominant matrix (~12)
  int expected_niters_prec; //!< expected number of test CG iterations with preconditioning and with diagonally dominant matrix (~1-2)
  int niters_max_no_prec; //!< maximum number of test CG iterations without predictitioner
  int niters_max_prec; //!< maximum number of test CG iterations without predictitioner
  Datatype normr; //!< residual norm achieved during test CG iterations
};

/*!
  Test the correctness of the Preconditined CG implementation by using a system matrix with a dominant diagonal.

  @param[in]    geom The description of the problem's geometry.
  @param[in]    A    The known system matrix
  @param[in]    data the data structure with all necessary CG vectors preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate solution
  @param[out]   testcg_data the data structure with the results of the test including pass/fail information

  @return Returns zero on success and a non-zero value otherwise.

  @see CG()
 */
template<class Datatype>
int TestCG(SparseMatrix<Datatype> & A, CGData<Datatype> & data, Vector<Datatype> & b, Vector<Datatype> & x, TestCGData<Datatype> & testcg_data) {


  // Use this array for collecting timing information
  std::vector< double > times(8,0.0);
  // Temporary storage for holding original diagonal and RHS
  Vector<Datatype> origDiagA, exaggeratedDiagA, origB;
  Geometry geom = *A.geom;
  InitializeVector(origDiagA, geom.nx, geom.ny, geom.nz);
  InitializeVector(exaggeratedDiagA, geom.nx, geom.ny, geom.nz);
  InitializeVector(origB, geom.nx, geom.ny, geom.nz);
  CopyMatrixDiagonal(A, origDiagA);
  CopyVector(origDiagA, exaggeratedDiagA);
  CopyVector(b, origB);

  // Modify the matrix diagonal to greatly exaggerate diagonal values.
  // CG should converge in about 10 iterations for this problem, regardless of problem size
  for (local_int_t i=0; i< A.localNumberOfRows; ++i) {
    global_int_t globalRowID = A.localToGlobalMap[i];
    if (globalRowID<9) {
      Datatype scale = (globalRowID+2)*1.0e6;
      ScaleVectorValue(exaggeratedDiagA, i, scale);
      ScaleVectorValue(b, i, scale);
    } else {
      ScaleVectorValue(exaggeratedDiagA, i, 1.0e6);
      ScaleVectorValue(b, i, 1.0e6);
    }
  }
  ReplaceMatrixDiagonal(A, exaggeratedDiagA);

  int niters = 0;
  Datatype normr = 0.0;
  Datatype normr0 = 0.0;
  int maxIters = 50;
  int numberOfCgCalls = 2;
  Datatype tolerance = 1.0e-12; // Set tolerance to reasonable value for grossly scaled diagonal terms
  testcg_data.expected_niters_no_prec = 12; // For the unpreconditioned CG call, we should take about 10 iterations, permit 12
  testcg_data.expected_niters_prec = 2;   // For the preconditioned case, we should take about 1 iteration, permit 2
  testcg_data.niters_max_no_prec = 0;
  testcg_data.niters_max_prec = 0;
  for (int k=0; k<2; ++k) { // This loop tests both unpreconditioned and preconditioned runs
    int expected_niters = testcg_data.expected_niters_no_prec;
    if (k==1) expected_niters = testcg_data.expected_niters_prec;
    for (int i=0; i< numberOfCgCalls; ++i) {
      ZeroVector(x); // Zero out x
      int ierr = CG(A, data, b, x, maxIters, tolerance, niters, normr, normr0, &times[0], k==1);
      if (ierr) HPCG_fout << "Error in call to CG: " << ierr << ".\n" << endl;
      if (niters <= expected_niters) {
        ++testcg_data.count_pass;
      } else {
        ++testcg_data.count_fail;
      }
      if (k==0 && niters>testcg_data.niters_max_no_prec) testcg_data.niters_max_no_prec = niters; // Keep track of largest iter count
      if (k==1 && niters>testcg_data.niters_max_prec) testcg_data.niters_max_prec = niters; // Same for preconditioned run
      if (A.geom->rank==0) {
        HPCG_fout << "Call [" << i << "] Number of Iterations [" << niters <<"] Scaled Residual [" << normr/normr0 << "]" << endl;
        if (niters > expected_niters)
          HPCG_fout << " Expected " << expected_niters << " iterations.  Performed " << niters << "." << endl;
      }
    }
  }

  // Restore matrix diagonal and RHS
  ReplaceMatrixDiagonal(A, origDiagA);
  CopyVector(origB, b);
  // Delete vectors
  DeleteVector(origDiagA);
  DeleteVector(exaggeratedDiagA);
  DeleteVector(origB);
  testcg_data.normr = normr;

  return 0;
}

#endif  // TESTCG_HPP
