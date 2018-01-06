% function OT problem with random data
% min sum(sum(C.*x))
clc; clear

% cost matrix
testResolution = 3;
costData = load('data/cost.mat');
C = costData.data{testResolution};

pic1Data = load('data/pic1.mat');
k = pic1Data.picN;
fprintf('Resolution %dx%d\n', pic1Data.resolutions(testResolution), pic1Data.resolutions(testResolution));

x0 = [];

errfun = @(x1, x2) norm(x1-x2)/(1+norm(x1));
draw = false;

%% test query
query = {...
%     {@dir_mosek, 'dir_mosek_simplex', struct('method', 'simplex')}; ...
%     {@transimplexWrapper, 'tran_simplex_normal', struct('method', 'normal')}; ...
    {@transimplexWrapper, 'tran_simplex_improved', struct('method', 'improved')}; ...
    {@transimplexWrapper, 'tran_simplex_shielding', struct('method', 'shielding')}; ...
%     {@transimplexWrapper, 'tran_simplex_shortlist', struct('method', 'shortlist')}; ...
    {@transimplexWrapper, 'tran_simplex_multiscale', struct('method', 'multiscale')}; ...
};
caseN = size(query, 1);

results = cell(k*(k-1)/2, caseN);
itr = 1;
for findex = 1:k
    mu = pic1Data.data{findex, testResolution};
    for bindex = findex+1:k
        nu = pic1Data.data{bindex, testResolution};
        for i = 1:caseN
            opts = query{i}{3};
            opts.draw = draw;
            tic;
            [results{itr, i}.x, results{itr, i}.out, results{itr, i}.val] = query{i}{1}(x0, C, mu, nu, opts);
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