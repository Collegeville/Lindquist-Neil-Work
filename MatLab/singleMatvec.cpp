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
	
	const float *x  = ((float*)mxGetData(prhs[3]))-1;

    const size_t Am = mxGetM(prhs[0])-1;
    
    mxArray *b = mxCreateUninitNumericMatrix(Am, 1, mxSINGLE_CLASS, mxREAL);
    plhs[0] = b;
    float *rawB = (float*)mxGetData(b);    
    
    //0-indexed row-index
    uint32_t row = Am;
    //0-indexed el-index of the value "before" i (i+1)
    uint32_t start = Ar[row]-1;
    while(row-- > 0){
        double temp = 0;
        //0-indexed el-index of i's last value
        uint32_t end = Ar[row]-1;
        //i is 1-indexed el-index
        for(uint32_t i = start; i-- > end;){
        	//i switches to 0-index el-index
            temp += (double)Av[i] * (double)x[Ac[i]];
        }
        rawB[row] = (float)temp;
        start = end;
    }
}


#ifdef __cplusplus
}
#endif
