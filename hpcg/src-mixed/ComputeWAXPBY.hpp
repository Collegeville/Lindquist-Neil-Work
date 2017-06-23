
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

#ifndef COMPUTEWAXPBY_HPP
#define COMPUTEWAXPBY_HPP
#include "Vector.hpp"
template<class datatype1, class datatype2, class datatype3>
int ComputeWAXPBY(const local_int_t n, const double alpha,
    const Vector<datatype1> & x, const double beta, const Vector<datatype2> & y,
    Vector<datatype3> & w, bool & isOptimized);
#endif // COMPUTEWAXPBY_HPP
