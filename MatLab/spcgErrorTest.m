function spcgErrorTest(continuation)

if nargin < 1
	continuation = true;
end


index = UFget();


%%%%%% single system %%%%%%%%
problemList = [1580, 1581, 1582, 1583, 1584, 1585, 1853, 1909, 1919, 2283];
problemAlphas=[ .05   .05   .01   .01   .01   .01   .01   .05    0    .01];

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
	if problemAlphas(i) == 0
		L1 = ichol(prob.A);
	else
		L1 = ichol(prob.A, struct('type', 'ict', 'droptol', 1e-3, 'diagcomp', problemAlphas(i)));
	end
	
	L1T = L1';

    b = prob.b;
	sb = single(full(b));

	[x, doubleFlags] = pcg(prob.A, b, 1e-6, 100, L1, L1T);
	sx = spcg(sA, sb, 1e-6, 100, L1, L1T);

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
L1 = ichol(prob.A, struct('type', 'ict', 'droptol', 1e-3, 'diagcomp', .05));
L1T = L1';


if continuation && exist('pcg results - 1850.mat', 'file') == 2
	load('pcg results - 1850.mat')
	disp(['loaded ' num2str(i) ' previous calculations']);
else
	doubleErr = zeros(problemCount, 1);
	singleErr = zeros(problemCount, 1);
	i = 1;
end
while i <= problemCount
    b = prob.b(:, 1);
    sb = single(full(b));
    
    [x, f] = pcg(prob.A, b, 1e-6, 100, L1, L1T);
    sx = spcg(sA, sb, 1e-06, 100, L1, L1T);
    
    doubleErr(i) = norm(b-prob.A*x);
    singleErr(i) = norm(b-prob.A*double(sx));
    
    disp(i);
    save('pcg results - 1850.mat', 'doubleErr', 'singleErr', 'i')

	i = i+1;
end