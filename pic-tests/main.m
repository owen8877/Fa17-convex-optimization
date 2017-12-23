% function OT problem with random data
% min sum(sum(C.*x))
clc; clear

% cost matrix
costData = load('data/cost.mat');
C = costData.data{3};

pic1Data = load('data/pic1.mat');
mu = pic1Data.data{1, 3};
nu = pic1Data.data{2, 3};
fprintf('Resolution %dx%d\n', pic1Data.resolutions(3), pic1Data.resolutions(3));

x0 = [];

errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
draw = false;

%% test query
query = {...
%     {@dir_mosek, 'dir_mosek_simplex', struct('method', 'simplex')}; ...
%     {@transimplexWrapper, 'tran_simplex_normal', struct('method', 'normal')}; ...
    {@transimplexWrapper, 'tran_simplex_improved', struct('method', 'improved')}; ...
    {@transimplexWrapper, 'tran_simplex_shortlist', struct('method', 'shortlist')}; ...
};
caseN = size(query, 1);

results = cell(1, caseN);
for i = 1:caseN
    opts = query{i}{3};
    opts.draw = draw;
    tic;
    [results{i}.x, results{i}.out, results{i}.val] = query{i}{1}(x0, C, mu, nu, opts);
    results{i}.t = toc;
end

%% print
cpu = cellfun(@(c) c.t, results)';
funcval = cellfun(@(c) full(c.out), results)';
error = cellfun(@(c) errfun(results{1}.x, c.x), results)';
method = cellfun(@(c) c{2}, query, 'UniformOutput', false)';
resTable = table(cpu, error, funcval, 'RowNames', method);

disp(resTable)

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