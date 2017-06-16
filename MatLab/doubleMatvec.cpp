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
    
    //0-indexed row-index
    mwIndex row = Am;
    
    mwIndex start = Ar[row];
    while(row-- > 0){
        double temp = 0;
        //0-indexed el-index of i's last value
        mwIndex end = Ar[row];
        //i is 1-indexed el-index
        for(mwIndex i = start; i-- > end;){
        	//i switches to 0-index el-index
            temp += Av[i] * x[Ac[i]];
        }
        rawB[row] = temp;
        start = end;
    }
}


#ifdef __cplusplus
}
#endif
