function matvecSpeedTest

A = rand(5000);
x = rand(5000, 1);
sA = single(A);
sx = single(x);

%single data, double arithmetic
mixedMult = @() matvecMixed(sA, sx);
%double data, auto arithmetic
doubleMult = @() matvecSame(A, x);
%double data, double arithmetic
doubleMult2 = @() matvecMixed(A, x);
%single data, auto arithmetic
singleMult = @() matvecSame(sA, sx);
%single data, single arithmetic
singleMult2 = @() matvecSingle(sA, sx);

eTime = timeit(doubleMult2)
sTime = timeit(singleMult)
dTime = timeit(doubleMult)
tTime = timeit(singleMult2)
mTime = timeit(mixedMult)