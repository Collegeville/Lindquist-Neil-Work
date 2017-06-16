/*
 * Computes the residual of b-A*x
 * r = residual(b, Ar, Ac, Av, x)
 * b and x should be single precision
 * Ar, Ac, Av should be the sparseSingle components
 * 
 * r is double precision
 */

#include <cstddef>
#include <cstdint>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]){

	if(nrhs != 5){
		mexErrMsgIdAndTxt("Neil:residual:InvalidArgCount", "Exactly 5 arguments are required");
	}

	const mxArray *b  = prhs[0];
	const float *rawB = (float*)mxGetData(b);

	//contains 1-indexed references
	const mxArray *Ar = prhs[1];
	const uint32_t *rawAr = (uint32_t*)mxGetData(Ar);

	//contains 1-indexed references
	const mxArray *Ac = prhs[2];
	const uint32_t *rawAc = (uint32_t*)mxGetData(Ac);

	const mxArray *Av = prhs[3];
	const float *rawAv = (float*)mxGetData(Av);

	const mxArray *x  = prhs[4];
	const float *rawX = ((float*)mxGetData(x))-1;

	const size_t m   = mxGetM(Ar)-1;

	mxArray *r = mxCreateUninitNumericMatrix(m, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[0] = r;
	double *rawR = mxGetPr(r);
    
    //0-indexed row-index
    uint32_t row = m;
    //0-indexed el-index of the value "before" i (i+1)
    uint32_t start = rawAr[row]-1;
    while(row-- > 0){
        double temp = 0;
        //0-indexed el-index of i's last value
        uint32_t end = rawAr[row]-1;
        //i is 1-indexed el-index
        for(uint32_t i = start; i-- > end;){
        	//i switches to 0-index el-index
            temp += (double)rawAv[i] * (double)rawX[rawAc[i]];
        }
        rawR[row] = (double)rawB[row] - temp;
        start = end;
    }
}

#ifdef __cplusplus
}
#endif