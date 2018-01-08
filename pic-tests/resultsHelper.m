function resultsHelper(R, range, method)
    mini = [];
    times = [];
    averagedTime = [];
    for testCase = range
        results = R{testCase};
        [tests, methodN] = size(results);
        
        out = cellfun(@(c) c.out, results);
        time = cellfun(@(c) c.t, results);
        mini = [mini; (out - min(out, [], 2)) ./ min(out, [], 2)];
        times = [times; time];
        averagedTime = [averagedTime; mean(time)];
    end
    figure
    boxplot(mini, 'Labels', method)
    ylabel('Relative error')
    set(gca, 'YScale', 'log')
    figure
    boxplot(times, 'Labels', method)
    ylabel('Run time (in seconds)')
end