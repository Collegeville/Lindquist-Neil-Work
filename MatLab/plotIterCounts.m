[sIter, dIter] = getIterCounts;
percentIncreases = 100*((sIter-dIter)./dIter);
avgPercentIncrease = mean(percentIncreases);


%reset the plot
clf('reset')
hold on

ylabel('Number of iterations')
plot(dIter)
plot(sIter)


yyaxis right
ylabel('percent increase')
plot(percentIncreases)
plot([1,length(dIter)],[avgPercentIncrease avgPercentIncrease])


xlabel('problem number')
legend('double precision iterations', 'mixed precision iterations', ...
    'percent increase', 'average percent increase')