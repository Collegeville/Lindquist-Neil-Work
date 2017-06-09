/*
 * Creates a transposed, double precision sparse array from
 * the component vectors of a sparseSingle, with vals cast to double
 *
 * sp = sparseCast_(m, n, rows, cols, vals);
 */

#include <cstdint>
#include <cstring>
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif
    
    
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    uint32_t m   = *(uint32_t*)mxGetData(prhs[0]);
    uint32_t n   = *(uint32_t*)mxGetData(prhs[1]);
    
    /*0 indexed, elements refer to 1 indexing*/
    const uint32_t *rowStarts = (uint32_t*)mxGetData(prhs[2]);
    /*0 indexed, elements refer to 1 indexing*/
    const uint32_t *colIndexes = (uint32_t*)mxGetData(prhs[3]);
    const float *vals = (float*)mxGetData(prhs[4]);
    
    mwSize nnz = mxGetM(prhs[4]);
    
    mxArray *sparse = mxCreateSparse(n, m, nnz, mxREAL);
    
    mwIndex *sparseRows = mxGetIr(sparse);
    mwIndex *sparseCols = mxGetJc(sparse);
    double  *sparseVals = mxGetPr(sparse);
    
    for(uint32_t i = 0; i<=m; i++){
        sparseCols[i] = (mwIndex)rowStarts[i]-1;
    }
    
    for(mwIndex i = 0; i<nnz; i++){
        sparseRows[i] = (mwIndex)colIndexes[i]-1;
        sparseVals[i] = (double)vals[i];
    }
    
    plhs[0] = sparse;
}
    
    
#ifdef __cplusplus
}
#endif
