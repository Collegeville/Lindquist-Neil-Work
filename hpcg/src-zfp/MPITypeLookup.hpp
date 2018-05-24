
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

#ifndef MPITYPELOOKUP_HPP
#define MPITYPELOOKUP_HPP

#include <mpi.h>

//credit https://stackoverflow.com/a/35869988/6353993

struct MPITypeWrapper {
    MPI_Datatype mpitype;
    int count;
};


/*!
  Gets the MPI datatype and the count for the given value.
  
  Unknown datatypes are given a type of MPI_Byte with a count equal to their size

  @param[in] t the value that will be sent

  @return returns MPITypeWrapper containing the type and count
*/
template<typename T> MPITypeWrapper getMPIType(T t) { MPITypeWrapper wrapper={ MPI_BYTE, sizeof(T) }; return wrapper; };
template<> MPITypeWrapper getMPIType(signed char t) { MPITypeWrapper wrapper={ MPI_SIGNED_CHAR, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(unsigned char t) { MPITypeWrapper wrapper={ MPI_UNSIGNED_CHAR, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(short t) { MPITypeWrapper wrapper={ MPI_SHORT, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(unsigned short t) { MPITypeWrapper wrapper={ MPI_UNSIGNED_SHORT, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(int t) { MPITypeWrapper wrapper={ MPI_INT, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(unsigned int t) { MPITypeWrapper wrapper={ MPI_UNSIGNED, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(long t) { MPITypeWrapper wrapper={ MPI_LONG, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(unsigned long t) { MPITypeWrapper wrapper={ MPI_UNSIGNED_LONG, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(float t) { MPITypeWrapper wrapper={ MPI_DOUBLE, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(double t) { MPITypeWrapper wrapper={ MPI_DOUBLE, 1 }; return wrapper; };
template<> MPITypeWrapper getMPIType(long double t) { MPITypeWrapper wrapper={ MPI_LONG_DOUBLE, 1 }; return wrapper; };

#endif 
