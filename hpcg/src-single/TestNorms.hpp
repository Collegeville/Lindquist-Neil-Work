
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
 @file TestNorms.hpp

 HPCG data structure
 */

#ifndef TESTNORMS_HPP
#define TESTNORMS_HPP


struct TestNormsData_STRUCT {
  float * values; //!< sample values
  float   mean;   //!< mean of all sampes
  float variance; //!< variance of mean
  int    samples;  //!< number of samples
  bool   pass;     //!< pass/fail indicator
};
typedef struct TestNormsData_STRUCT TestNormsData;

extern int TestNorms(TestNormsData & testnorms_data);

#endif  // TESTNORMS_HPP
