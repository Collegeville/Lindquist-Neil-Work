function [singleIters, doubleIters] = getIterCounts()

probIds = [1580, 1581, 1582, 1583, 1584, 1585, 1853, 1909, 1919, 2283];

probCount = length(probIds);

singleIters = zeros(probCount, 1);
doubleIters = zeros(probCount, 1);

for p = 1:probCount
    load(['pcg results - ' num2str(probIds(p)) '.mat'])
    singleIters(p) = sIter;
    doubleIters(p) = dIter;
end


i = 0;
load('pcg results - 1850.mat')
doubleIters(probCount+1:probCount+i-1) = sIter(1:i-1);
singleIters(probCount+1:probCount+i-1) = dIter(1:i-1);
