function [avgY, avg2] = errortest(iter)

if nargin < 1
    iter = 10000;
end

mat = [4, -1, -1, -1, -1];
ssp = sparseSingle(mat);
smat = single(mat);

avgY = zeros(10, 1);
avgErrS = zeros(10, 1);
avgErrM = zeros(10, 1);

for cntr = 1:iter
    r = rand(5, 1);
    
    for i = 1:10
        x = 1+r/(10^i);
        sx = single(x);
                
        d = mat*x;
        s = smat*sx;
        m = ssp*sx;
        
        errS = norm(d-double(s));
        errM = norm(d-double(m));
        
        avgY(i) = avgY(i) + log2(errS/errM);
        avgErrS(i) = avgErrS(i) + errS;
        avgErrM(i) = avgErrM(i) + errM;
    end
end
avgY = 2.^(avgY/iter);
avg2 = avgErrS./avgErrM;

plot(1:10, avgY, 1:10, avg2);