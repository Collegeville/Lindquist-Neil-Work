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
	const float *rawX = (float*)mxGetData(x);

	const size_t m   = mxGetM(Ar)-1;
	const size_t n   = mxGetM(x);
	const size_t nnz = mxGetM(Av);

	mxArray *r = mxCreateUninitNumericMatrix(m, 1, mxDOUBLE_CLASS, mxREAL);
	double *rawR = mxGetPr(r);

	//0-indexed value indexes
	uint32_t nextEnd = rawAr[1]-1;
	//0-indexed row indexes
	uint32_t row = 0;
	double temp;
	//0-indexed value indexes
	uint32_t i = 0;
	while(i < nnz){
		if(i == nextEnd){
			rawR[row] = rawB[row]-temp;
			temp = 0;
			row++;
			nextEnd = rawAr[row+1]-1;
		}
		temp += (double)rawAv[i]*(double)rawX[rawAc[i]-1];
		i++;
	}
	rawR[row] = rawB[row]-temp;

	plhs[0] = r;
}

#ifdef __cplusplus
}
#endif