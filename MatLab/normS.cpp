/* 
 * Computes the Euclidean Norm of two single precision floats
 * and returns the result in double precision.
 * 
 * n = normS(a) should be equivalent to n = norm(double(a));
 */

#include <cmath>
#include <cstddef>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[]){

	if(nrhs != 1){
		mexErrMsgIdAndTxt("Neil:normS:InvalidArgCount", "Exactly 1 arguments are required");
	}

    size_t length = mxGetM(prhs[0]);
    
    float *rawA = (float*)mxGetData(prhs[0]);
    
    double sum = 0;
    for(size_t i = 0; i<length; i++){
        double entry = (double)rawA[i];
        sum += entry*entry;
    }
    
    plhs[0] = mxCreateDoubleScalar(sqrt(sum));
}

#ifdef __cplusplus
}
#endif