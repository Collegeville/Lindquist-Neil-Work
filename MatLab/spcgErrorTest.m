
index = UFget();

prob = UFget(1850, index);
disp([prob.name ' - ' num2str(prob.id)]);

problemCount = size(prob.b, 2);

%40543 b's
sA = sparseSingle(prob.A);
L1 = ichol(prob.A, struct('type', 'ict', 'droptol', 1e-3, 'diagcomp', .05));
L1T = L1';

doubleErr = zeros(problemCount, 1);
singleErr = zeros(problemCount, 1);
disp('setup complete');
for i = 1:problemCount
    b = prob.b(:, 1);
    sb = single(full(b));
    
    [x, f] = pcg(prob.A, b, 1e-6, 100, L1, L1T);
    sx = spcg(sA, sb, 1e-06, 100, L1, L1T);
    
    doubleErr(i) = norm(b-prob.A*x);
    singleErr(i) = norm(b-prob.A*double(sx));
    
    disp(i);
    save('pcg results - 1850.m', 'doubleErr', 'singleErr', 'i')
end