
[sErr, dErr, nB] = getErrors;


%reset the plot
clf('reset')


semilogy(dErr./nB)
hold on
semilogy(sErr./nB)
semilogy([1,length(dErr)],[1e-6 1e-6], 'k-')
semilogy([1,length(dErr)],[1e-5 1e-5], 'k-')
%semilogy([10.5,10.5],[1e-4 1e-7], 'k-')


ylabel('Relative Error of Result')
xlabel('Problem Number')
legend('Double Precision Error', 'Mixed Precision Error')