/*
 * Computes the dot product between two single precision floats
 * d = dotProdSS(a, b)
 * Is equivalent to d = double(a')*double(b)
 */

#include <cstddef>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 2){
        mexErrMsgIdAndTxt("Neil:dotProdSS:InvalidArgCount", "Exactly 2 arguments are required");
    }
    
    size_t length = mxGetM(prhs[0]);
    if(mxGetM(prhs[1]) != length){
        mexErrMsgIdAndTxt("Neil:dotProdSS:MismatchedLengths", "Vectors must be the same size");
    }
    
    float *rawA = (float*)mxGetData(prhs[0]);
    float *rawB = (float*)mxGetData(prhs[1]);
    
    double sum = 0;
    for(size_t i = 0; i<length; i++){
        sum += ((double)rawA[i])*((double)rawB[i]);
    }
    
    plhs[0] = mxCreateDoubleScalar(sum);
}


#ifdef __cplusplus
}
#endif