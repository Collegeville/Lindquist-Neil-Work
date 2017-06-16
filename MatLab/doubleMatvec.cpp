/*
 * Computes the product of a the transpose of a sparse matrix and dense vector
 * in double precision.
 * 
 * b = singleMatvec(AT, x)
 */


#include <cstdint>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
	const mxArray *A = prhs[0];
    
    //0-indexed
    const mwIndex *Ac = mxGetIr(A);
    //0-indexed
    const mwIndex *Ar = mxGetJc(A);
    const double  *Av = mxGetPr(A);
	
	//0 indexed
	const double *x  = mxGetPr(prhs[1]);
    
    const size_t Am = mxGetN(A);
    
    mxArray *b = mxCreateUninitNumericMatrix(Am, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[0] = b;
    double *rawB = mxGetPr(b);    
    
    double temp = 0;
    //0-indexed row-index
    mwIndex row = Am-1;
    //0-indexed el-index
    mwIndex nextEnd = Ar[row];
    //1-indexed el-index
    mwIndex i = Ar[row+1];
    while(i > 0){
        if(i <= nextEnd){
            rawB[row] = temp;
            --row;
            nextEnd = Ar[row];
            temp = 0;
        }else{
            --i;
            temp += Av[i] * x[Ac[i]];
        }
    }
    rawB[row] = temp;
    while(row > 0){
        rawB[--row] = 0;
    }
}


#ifdef __cplusplus
}
#endif
