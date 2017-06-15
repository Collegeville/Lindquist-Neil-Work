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
	
	const float *x  = (float*)mxGetData(prhs[3])-1;

    const size_t Am = mxGetM(prhs[0])-1;
    
    mxArray *b = mxCreateUninitNumericMatrix(Am, 1, mxSINGLE_CLASS, mxREAL);
    plhs[0] = b;
    float *rawB = (float*)mxGetData(b);    
    
    double temp = 0;
    //0-indexed row-index
    uint32_t row = Am-1;
    //1-indexed el-index
    uint32_t nextEnd = Ar[row];
    //1-indexed el-index
    uint32_t i = Ar[row+1]-1;
    while(i > 0){
        if(i < nextEnd){
            rawB[row] = (float)temp;
            --row;
            nextEnd = Ar[row];
            temp = 0;
        }else{
            --i;
            temp += (double)Av[i] * (double)x[Ac[i]];
        }
    }
    rawB[row] = (float)temp;
    while(row > 0){
        rawB[--row] = 0;
    }
}


#ifdef __cplusplus
}
#endif