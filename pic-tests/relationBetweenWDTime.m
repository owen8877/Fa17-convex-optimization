function resultsHelper(R, range, method)
    mini = [];
    times = [];
    averagedTime = [];
    wd = cell2mat(cellfun(@(c) mean(mean(cellfun(@(d) d.out, c))), R, 'UniformOutput', false));
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
    hold on
    for i = 1:size(averagedTime, 2)
        scatter(averagedTime(:, i), wd);
    end
    legend(method);
end