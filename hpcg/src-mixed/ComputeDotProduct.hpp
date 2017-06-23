
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
template<class datatype1, class datatype2>
int ComputeDotProduct(const local_int_t n, const Vector<datatype1> & x,
    const Vector<datatype2> & y, double & result, double & time_allreduce,
    bool & isOptimized);

#endif // COMPUTEDOTPRODUCT_HPP
