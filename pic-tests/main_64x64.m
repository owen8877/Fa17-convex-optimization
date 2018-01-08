% function OT problem with random data
% min sum(sum(C.*x))
% clc; clear

% cost matrix
testResolution = 5;
costData = load('data/cost.mat');
C = costData.data{testResolution};
errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
draw = false;
R = cell(10, 1);

%% test query
testRange = 1:1;
query = {...
%     {@dir_mosek, 'mosek', struct('method', 'simplex')}; ...
%     {@dir_cplex, 'cplex_lp', []}; ...
    {@transimplexWrapper, 'cplex_net', struct('method', 'cplex')}; ...
%     {@transimplexWrapper, 'normal', struct('method', 'normal')}; ...
    {@transimplexWrapper, 'improved', struct('method', 'improved')}; ...
    {@transimplexWrapper, 'shielding', struct('method', 'shielding')}; ...
%     {@transimplexWrapper, 'shortlist', struct('method', 'shortlist')}; ...
    {@transimplexWrapper, 'multiscale', struct('method', 'multiscale')}; ...
};
caseN = size(query, 1);

for testCase = testRange
    picData = load(['data/pic', int2str(testCase), '.mat']);
    k_ = picData.picN;
    k = min(3, k_); % real test pic numbers
    fprintf('Resolution %dx%d\n', picData.resolutions(testResolution), picData.resolutions(testResolution));

    results = cell(k*(k-1)/2, caseN);
    itr = 1;
    for findex = 1:k
        mu = picData.data{findex, testResolution};
        for bindex = findex+1:k
            nu = picData.data{bindex, testResolution};
            for i = 1:caseN
                opts = query{i}{3};
                opts.draw = draw;
                tic;
                [results{itr, i}.x, results{itr, i}.out, results{itr, i}.val] = query{i}{1}([], C, mu, nu, opts);
                results{itr, i}.t = toc;
            end
            itr = itr + 1;
        end
    end

    %% print
    cpu = cellfun(@(c) c.t, results)';
    cpu = cpu(:, 1);
    funcval = cellfun(@(c) full(c.out), results)';
    funcval = funcval(:, 1);
    error = cellfun(@(c) errfun(results{1}.x, c.x), results)';
    error = error(:, 1);
    method = cellfun(@(c) c{2}, query, 'UniformOutput', false)';
    resTable = table(cpu, error, funcval, 'RowNames', method);
    fprintf('The first test case:\n')
    disp(resTable)

    cpu = cellfun(@(c) c.t, results)';
    cpu = cpu(:, 1);
    method = cellfun(@(c) c{2}, query, 'UniformOutput', false)';
    resTable = table(cpu, 'RowNames', method);
    fprintf('Average test case:\n')
    disp(resTable)
    R{testCase} = results;
end
% resultsHelper(R, 1:1, method)

% if draw
%     globalMin = min([results{2}.val results{3}.val results{4}.val results{5}.val]);
%     semilogy(1:numel(results{2}.val), results{2}.val-globalMin)
%     hold on
%     semilogy(1:numel(results{3}.val), results{3}.val-globalMin)
%     semilogy(1:numel(results{4}.val), results{4}.val-globalMin)
%     semilogy(1:numel(results{5}.val), results{5}.val-globalMin)
%     title('Four opt methods')
%     legend({'momentum conti neste', 'adagrad conti neste', 'rmsprop conti', 'adam conti'})
%     xlabel('k')
%     ylabel('f(x_k)-f(x^*)')
%     set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.5);
% end