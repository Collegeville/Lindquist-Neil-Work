/*
 * Computes the product of a sparse matrix and dense vector
 * in single precision.
 * 
 * b = singleMatvec(Ar, Ac, Av, x)
 */


#include <cstdint>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
	const uint32_t *Ar = (uint32_t*)mxGetData(prhs[0]);
	const uint32_t *Ac = (uint32_t*)mxGetData(prhs[1]);
	const float *Av = (float*)mxGetData(prhs[2]);
	
	const float *x  = (float*)mxGetData(prhs[3]);

    const size_t Am   = mxGetM(prhs[0])-1;
    
    mxArray *b = mxCreateUninitNumericMatrix(Am, 1, mxSINGLE_CLASS, mxREAL);
    plhs[0] = b;
    float *rawB = (float*)mxGetData(b);    
    
    uint32_t nextEnd = Ar[1]-1;
    double temp = 0;
    uint32_t i = 0;
    uint32_t row = 0;
    while(row != Am){
        if(i == nextEnd){
            rawB[row] = (float)temp;
            ++row;
            nextEnd = Ar[row+1]-1;
            temp = 0;
        }else{
            temp += (double)Av[i] * (double)x[Ac[i]-1];
            ++i;
        }
    }
    rawB[row] = (float)temp;
}


#ifdef __cplusplus
}
#endif