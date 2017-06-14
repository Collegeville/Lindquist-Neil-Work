/*
 * Computes the product of a sparse matrix and dense vector
 * in single precision.
 * 
 * b = singleMatvec(Am, An, Annz, Ar, Ac, Av, x)
 */

#include <cstdint>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

	uint32_t Am   = *(uint32_t*)mxGetData(prhs[0]);
	uint32_t An   = *(uint32_t*)mxGetData(prhs[1]);
	uint32_t Annz = *(uint32_t*)mxGetData(prhs[2]);

	const uint32_t *Ar = (uint32_t*)mxGetData(prhs[3]);
	const uint32_t *Ac = (uint32_t*)mxGetData(prhs[4]);
	const float *Av = (float*)mxGetData(prhs[5]);
	
	const float *x  = (float*)mxGetData(prhs[6]);

    float *rawB = (float*)mxCalloc(Am, sizeof(float));

    
    uint32_t nextEnd = Ar[1]-1;
    uint32_t row = 0;
    double temp = 0;
    uint32_t i = 0;
    while(i < Annz){
        if(i == nextEnd){
            rawB[row] = (float)temp;
            temp = 0;
            row++;
            nextEnd = Ar[row+1]-1;
        }
        temp += (double)Av[i] * (double)x[Ac[i]-1];
        i++;
    }
    rawB[row] = (float)temp;
    
    mxArray *b = mxCreateUninitNumericMatrix(Am, 1, mxSINGLE_CLASS, mxREAL);
    mxFree(mxGetData(b));
    mxSetData(b, rawB);
    
    plhs[0] = b;
}


#ifdef __cplusplus
}
#endif