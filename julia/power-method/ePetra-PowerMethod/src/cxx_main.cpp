//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#define EPETRA_MPI


#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <cmath>
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Version.h"

// prototype
int power_method(Epetra_CrsMatrix& A, double & lambda, int niters, double tolerance,
		 bool verbose);


int main(int argc, char *argv[])
{
  int ierr = 0, i;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);

  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  Epetra_SerialComm Comm;

#endif

  int MyPID = Comm.MyPID();
  bool verbose = (MyPID==0);

  if (verbose){
    std::cout << Epetra_Version() << std::endl << std::endl;
    std::cout << Comm << std::endl;
  }


  Epetra_CrsMatrix * Aptr;
  EPETRA_CHK_ERR(EpetraExt::MatrixMarketFileToCrsMatrix64("power-method/matrix.mm", Comm, Aptr, 0, verbose));

  Epetra_CrsMatrix & A = *Aptr;

  // variable needed for iteration
  double lambda = 0.0;
  int niters = (int) A.NumGlobalRows64()*10;

  double tolerance = 1.0e-2;


  // Iterate
  Epetra_Time timer(Comm);
  ierr += power_method(A, lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();


  if (verbose){
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "Total Time for first (warm up) solve = " << elapsed_time << std::endl<< std::endl;
  }


  // Iterate (again)
  lambda = 0.0;
  timer.ResetStartTime();
  ierr += power_method(A, lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();

  if (verbose){
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "Total time for first (warmed) solve = " << elapsed_time << std::endl<< std::endl;
  }

  // Increase diagonal dominance
  if (verbose)
    std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n\n"
              << std::endl;

  if (A.MyGlobalRow(0)) {
    int numvals = A.NumGlobalEntries(0);
    std::vector<double> Rowvals(numvals);
    std::vector<long long> Rowinds(numvals);
    A.ExtractGlobalRowCopy(0, numvals, numvals, &Rowvals[0], &Rowinds[0]); // Get A[0,0]
    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;

    A.ReplaceGlobalValues(0, numvals, &Rowvals[0], &Rowinds[0]);
  }

  // Iterate (again)
  lambda = 0.0;
  timer.ResetStartTime();
  ierr += power_method(A, lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();

  if (verbose){
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "Total time for second solve = " << elapsed_time << std::endl<< std::endl;
    std::cout << "\nErrors: " << ierr << std::endl << std::endl;
  }


  // Release all objects
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(Epetra_CrsMatrix& A, double &lambda, int niters,
        double tolerance, bool verbose) {

  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());

  // Fill z with random Numbers
  z.Random();

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;

  int iter = 0;
  for (; iter < niters; iter++)
  {
    z.Norm2(&normz); // Compute 2-norm of z
    q.Scale(1.0/normz, z);
    A.Multiply(false, q, z); // Compute z = A*q
    q.Dot(z, &lambda); // Approximate maximum eigenvalue
    if (iter%100==0 || iter+1==niters)
    {
      resid.Update(1.0, z, -lambda, q, 0.0); // Compute A*q - lambda*q
      resid.Norm2(&residual);
    }
    if (residual < tolerance) {
      ierr = 0;
      break;
    }
  }
  if (verbose) {
    std::cout << "Number of iterations: " << iter << std::endl;
  }
  return(ierr);
}
