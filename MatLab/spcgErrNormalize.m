function [normedErrList, normedErr1850] = spcgErrNormalize(process1850)


%don't process 1850 without explicit request
%it may be being processed in spcgErrorTest
if nargin < 1
    process1850 = false;
end

probIds = [1580, 1581, 1582, 1583, 1584, 1585, 1853, 1909, 1919, 2283];

normedErrList = zeros(length(probIds), 1);

for p = 1:length(probIds)
    load(['pcg results - ' num2str(probIds(p)) '.mat'])
    normedErr = singleErr/doubleErr;
    normedErrList(p) = normedErr;
    save(['pcg results - ' num2str(probIds(p)) '.mat'], 'doubleErr', 'singleErr', 'doubleFlags', 'normedErr')
    clear doubleErr singleErr doubleFlags normedErr
end

% disp(normedErrList)


if process1850
    i = 0;
    load('pcg results - 1850.mat')
    doubleErrTrunc = doubleErr(1:i-1);
    singleErrTrunc = singleErr(1:i-1);
    normedErr1850 = singleErrTrunc./doubleErrTrunc;
    save('pcg results - 1850.mat', 'doubleErr', 'singleErr', 'normedErr1850', 'i')
else
    normedErr1850 = zeros(0, 1);
end