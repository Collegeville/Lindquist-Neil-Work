
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

#ifndef EXCHANGEHALO_HPP
#define EXCHANGEHALO_HPP
// Compile this routine only if running with MPI
#ifndef HPCG_NO_MPI

#include <mpi.h>
#include "SparseMatrix.hpp"
#include "Vector.hpp"
#include "Geometry.hpp"
#include "ExchangeHalo.hpp"
#include "MPITypeLookup.hpp"
#include <cstdlib>

/*!
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated; on exit: the vector with non-local entries updated by other processors
 */
template<class Datatype>
void ExchangeHalo(const SparseMatrix<Datatype> & A, Vector<Datatype> & x) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  int num_neighbors = A.numberOfSendNeighbors;
  local_int_t * receiveLength = A.receiveLength;
  local_int_t * sendLength = A.sendLength;
  int * neighbors = A.neighbors;
  Datatype * sendBuffer = A.sendBuffer;
  local_int_t totalToBeSent = A.totalToBeSent;
  local_int_t * elementsToSend = A.elementsToSend;

  zfp::array3<Datatype>& xv = x.values;

  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //
  //  first post receives, these are immediate receives
  //  Do not wait for result to come, will do that at the
  //  wait call below.
  //

  int MPI_MY_TAG = 99;

  MPI_Request * request = new MPI_Request[num_neighbors];

  Datatype* x_external = new Datatype[localNumberOfRows];

  MPITypeWrapper mpitype = getMPIType(xv[0]);

  // Post receives first
  // TODO: Thread this loop
  for (int i = 0; i < num_neighbors; i++) {
    local_int_t n_recv = receiveLength[i]*mpitype.count;
    MPI_Irecv(x_external, n_recv, mpitype.mpitype, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
    x_external += n_recv;
  }


  //
  // Fill up send buffer
  //

  // TODO: Thread this loop
  for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xv[elementsToSend[i]];

  //
  // Send to each neighbor
  //

  // TODO: Thread this loop
  for (int i = 0; i < num_neighbors; i++) {
    local_int_t n_send = sendLength[i]*mpitype.count;
    MPI_Send(sendBuffer, n_send, mpitype.mpitype, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
    sendBuffer += n_send;
  }

  //
  // Complete the reads issued above
  //

  MPI_Status status;
  // TODO: Thread this loop
  for (int i = 0; i < num_neighbors; i++) {
    if ( MPI_Wait(request+i, &status) ) {
      std::exit(-1); // TODO: have better error exit
    }
  }

  //
  // Externals are at end of locals
  //
  for (int i = 0; i < localNumberOfRows; i++) {
    xv[i+localNumberOfRows] = x_external[i];
  }
  delete [] x_external;

  delete [] request;

  return;
}
#endif
#endif // EXCHANGEHALO_HPP
