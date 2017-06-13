
[sErr, dErr, nB] = getErrors;

hold off
semilogy(sErr./nB)
hold on
semilogy(dErr./nB)
semilogy([1,length(dErr)],[1e-6 1e-6], 'k-')
semilogy([1,length(dErr)],[1e-5 1e-5], 'k-')
semilogy([10.5,10.5],[1e-4 1e-7], 'k-')