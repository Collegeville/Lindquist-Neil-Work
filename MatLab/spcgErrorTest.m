function [errRatio, avgErr] = spcgErrorTest(iter, sz)
    
    if nargin < 1
        iter = 10;
    end
    if nargin < 2
        sz = 100;
    end
    
    errRatio = 0;
    avgErr = 0;

    for i = 1:iter
        A = sprandsym(sz, 3/sz, .1, 1);
        x = rand(sz, 1);
        b = A*x;
        
        As = sparseSingle(A);
        bs = single(b);
        
        xm = spcg(As, bs);
        xd = spcg(A, b);
        
        err = norm(double(xm)-x);
        errRatio = errRatio + log2(err/norm(xd-x));
        avgErr = avgErr + err;

    end
    
    errRatio = 2^(errRatio/iter);
    avgErr = avgErr/iter;