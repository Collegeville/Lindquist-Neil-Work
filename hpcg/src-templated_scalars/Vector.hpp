
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
 @file Vector.hpp

 HPCG data structures for dense vectors
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP
#include <cassert>
#include <cstdlib>
#include "Geometry.hpp"

template<class Datatype>
struct Vector {
  local_int_t localLength;  //!< length of local portion of the vector
  Datatype * values;          //!< array of values
  /*!
   This is for storing optimized data structures created in OptimizeProblem and
   used inside optimized ComputeSPMV().
   */
  Datatype * optimizationData;

};


/*!
  Initializes input vector.

  @param[in] v
  @param[in] localLength Length of local portion of input vector
 */
template<class Datatype>
inline void InitializeVector(Vector<Datatype> & v, local_int_t localLength) {
  v.localLength = localLength;
  v.values = new Datatype[localLength];
  v.optimizationData = 0;
  return;
}

/*!
  Fill the input vector with zero values.

  @param[inout] v - On entrance v is initialized, on exit all its values are zero.
 */
template<class Datatype>
inline void ZeroVector(Vector<Datatype> & v) {
  local_int_t localLength = v.localLength;
  Datatype * vv = v.values;
  for (int i=0; i<localLength; ++i) vv[i] = 0.0;
  return;
}
/*!
  Multiply (scale) a specific vector entry by a given value.

  @param[inout] v Vector to be modified
  @param[in] index Local index of entry to scale
  @param[in] value Value to scale by
 */
template<class Datatype1, class Datatype2>
inline void ScaleVectorValue(Vector<Datatype1> & v, local_int_t index, Datatype2 value) {
  assert(index>=0 && index < v.localLength);
  Datatype1 * vv = v.values;
  vv[index] *= value;
  return;
}
/*!
  Fill the input vector with pseudo-random values.

  @param[in] v
 */
template<class Datatype>
inline void FillRandomVector(Vector<Datatype> & v) {
  local_int_t localLength = v.localLength;
  Datatype * vv = v.values;
  for (int i=0; i<localLength; ++i) vv[i] = rand() / (Datatype)(RAND_MAX) + 1.0;
  return;
}
/*!
  Copy input vector to output vector.

  @param[in] v Input vector
  @param[in] w Output vector
 */
template<class Datatype1, class Datatype2>
inline void CopyVector(const Vector<Datatype1> & v, Vector<Datatype2> & w) {
  local_int_t localLength = v.localLength;
  assert(w.localLength >= localLength);
  Datatype1 * vv = v.values;
  Datatype2 * wv = w.values;
  for (int i=0; i<localLength; ++i) wv[i] = (Datatype2)vv[i];
  return;
}


/*!
  Deallocates the members of the data structure of the known system matrix provided they are not 0.

  @param[in] A the known system matrix
 */
template<class Datatype>
inline void DeleteVector(Vector<Datatype> & v) {

  delete [] v.values;
  v.localLength = 0;
  return;
}

#endif // VECTOR_HPP