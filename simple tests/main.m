% function OT problem with random data
% min sum(sum(C.*x))
clc; clear

%% generate data
n = 64;
m = 64;

% random cost
C = rand(m,n);
% distance cost

mu = rand(m, 1);
mu = mu / sum(mu);
nu = rand(n, 1);
nu = nu / sum(nu);
xSample = rand(m, n);
nu = sum(xSample, 1)';
nu = nu / sum(nu);
mu = sum(xSample, 2);
mu = mu / sum(mu);

x0 = rand(n*m, 1);

errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
draw = false;

%% test query
query = {...
    {@dir_mosek, 'dir_mosek_interior', struct('method', 'interior')}; ...
    {@dir_mosek, 'dir_mosek_simplex', struct('method', 'simplex')}; ...
%     {@admm, 'admm', struct('nesterov', true, 'tor', 1e-4)}; ...
%     {@admm, 'admm', struct('nesterov', false, 'tor', 5e-5)}; ...
    {@transimplex, 'transimplex', struct()}
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
funcval = cellfun(@(c) c.out, results)';
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