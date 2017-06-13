function [singleErrs, doubleErrs, normB] = getErrors()

index = UFget();

probIds = [1580, 1581, 1582, 1583, 1584, 1585, 1853, 1909, 1919, 2283];

probCount = length(probIds);

singleErrs = zeros(probCount, 1);
doubleErrs = zeros(probCount, 1);
normB = zeros(probCount, 1);

for p = 1:probCount
    load(['pcg results - ' num2str(probIds(p)) '.mat'])
    singleErrs(p) = singleErr;
    doubleErrs(p) = doubleErr;
    
    prob = UFget(probIds(p), index);
    normB(p) = norm(prob.b);
end


i = 0;
load('pcg results - 1850.mat')
doubleErrs(probCount+1:probCount+i-1) = doubleErr(1:i-1);
singleErrs(probCount+1:probCount+i-1) = singleErr(1:i-1);

prob = UFget(1850, index);
b = prob.b(:, 1:i-1);
normB(probCount+1:probCount+i-1) = sqrt(sum(b.^2,1));
