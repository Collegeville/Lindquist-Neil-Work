
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
 @file CGData.hpp

 HPCG data structure
 */

#ifndef CGDATA_HPP
#define CGDATA_HPP

#include "SparseMatrix.hpp"
#include "Vector.hpp"

template<class Datatype>
struct CGData {
  Vector<Datatype> r; //!< pointer to residual vector
  Vector<Datatype> z; //!< pointer to preconditioned residual vector
  Vector<Datatype> p; //!< pointer to direction vector
  Vector<Datatype> Ap; //!< pointer to Krylov vector
};

/*!
 Constructor for the data structure of CG vectors.

 @param[in]  A    the data structure that describes the problem matrix and its structure
 @param[out] data the data structure for CG vectors that will be allocated to get it ready for use in CG iterations
 */
template<class DatatypeA, class DatatypeData>
inline void InitializeSparseCGData(SparseMatrix<DatatypeA> & A, CGData<DatatypeData> & data) {
  Geometry geom = *A.geom;
  InitializeVector(data.r, geom.nx, geom.ny, geom.nz);
  InitializeVector(data.z, geom.nx, geom.ny, geom.nz);
  InitializeVector(data.p, geom.nx, geom.ny, geom.nz);
  InitializeVector(data.Ap, geom.nx, geom.ny, geom.nz);
  return;
}

/*!
 Destructor for the CG vectors data.

 @param[inout] data the CG vectors data structure whose storage is deallocated
 */
template<class Datatype>
inline void DeleteCGData(CGData<Datatype> & data) {

  DeleteVector (data.r);
  DeleteVector (data.z);
  DeleteVector (data.p);
  DeleteVector (data.Ap);
  return;
}

#endif // CGDATA_HPP

