
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

#include <cmath>

template<class Datatype>
struct TestNormsData {
  Datatype * values; //!< sample values
  Datatype   mean;   //!< mean of all sampes
  Datatype variance; //!< variance of mean
  int    samples;  //!< number of samples
  bool   pass;     //!< pass/fail indicator
};

/*!
  Computes the mean and standard deviation of the array of norm results.

  @param[in] testnorms_data data structure with the results of norm test

  @return Returns 0 upon success or non-zero otherwise
*/
template<class Datatype>
int TestNorms(TestNormsData<Datatype> & testnorms_data) {
 Datatype mean_delta = 0.0;
 for (int i= 0; i<testnorms_data.samples; ++i) mean_delta += (testnorms_data.values[i] - testnorms_data.values[0]);
 Datatype mean = testnorms_data.values[0] + mean_delta/(Datatype)testnorms_data.samples;
 testnorms_data.mean = mean;

 // Compute variance
 Datatype sumdiff = 0.0;
 for (int i= 0; i<testnorms_data.samples; ++i) sumdiff += (testnorms_data.values[i] - mean) * (testnorms_data.values[i] - mean);
 testnorms_data.variance = sumdiff/(Datatype)testnorms_data.samples;

 // Determine if variation is sufficiently small to declare success
 testnorms_data.pass = (testnorms_data.variance<1.0e-6);

 return 0;
}

#endif  // TESTNORMS_HPP
