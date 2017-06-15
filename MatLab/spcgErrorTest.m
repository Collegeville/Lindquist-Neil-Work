function spcgErrorTest(continuation)

if nargin < 1
	continuation = true;
end


index = UFget();


%%%%%% single system %%%%%%%%
problemList = [1580,   1581,  1582,  1583,  1584,  1585,  1853,  1909,  1919,  2283];
problemAlphas=[.0001   .0001  .0001  .0001  .0001  .0001  .0001  .0001  .0001   .01];
droptols     =[2.5e-5 2.5e-5 2.5e-5 2.5e-5 2.5e-5 2.5e-5 2.5e-5 1.5e-5  1e-3   1e-3];

%each has a problemCount of 1

for i = 1:length(problemList)
    
    p = problemList(i);
    
    if continuation
        if exist(['pcg results - ' num2str(p) '.mat'], 'file') == 2
            disp(['Problem ' num2str(p) ' already computed']);
            continue
        end
    end
    
	prob = UFget(p, index);

	disp([prob.name ' - ' num2str(prob.id)]);

	sA = sparseSingle(prob.A);
	L1 = ichol(prob.A, struct('type', 'ict', 'droptol', droptols(i), 'diagcomp', problemAlphas(i)));
	
	L1T = L1';

    b = prob.b;
	sb = single(full(b));

	[x, doubleFlags] = pcg(prob.A, b, 1e-6, 260, L1, L1T);
	sx = spcg(sA, sb, 1e-6, 260, L1, L1T);

	doubleErr = norm(b-prob.A*x);
    singleErr = norm(b-prob.A*double(sx));

	disp(p);
	save(['pcg results - ' num2str(p) '.mat'], 'doubleErr', 'singleErr', 'doubleFlags');	
end



%%%%%%%%%%% 1850 %%%%%%%%%%%%%
prob = UFget(1850, index);
disp([prob.name ' - ' num2str(prob.id)]);

problemCount = size(prob.b, 2);

%40543 b's
sA = sparseSingle(prob.A);
L1 = ichol(prob.A, struct('type', 'ict', 'droptol', 2e-5, 'diagcomp', .0001));
L1T = L1';


if continuation && exist('pcg results - 1850.mat', 'file') == 2
	load('pcg results - 1850.mat')
	disp(['loaded ' num2str(i-1) ' previous calculations']);
else
	doubleErr = zeros(problemCount, 1);
	singleErr = zeros(problemCount, 1);
	i = 1;
end
while i <= problemCount
    b = prob.b(:, i);
    sb = single(full(b));
    
    [x, ~] = pcg(prob.A, b, 1e-6, 400, L1, L1T);
    sx = spcg(sA, sb, 1e-06, 400, L1, L1T);
    
    doubleErr(i) = norm(b-prob.A*x);
    singleErr(i) = norm(b-prob.A*double(sx));
    
    disp(i);
    i = i+1;
    save('pcg results - 1850.mat', 'doubleErr', 'singleErr', 'i')
end