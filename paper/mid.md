## 凸优化（2017秋）选题报告

__陈子恒 1500010632__

本人选题为 第二题 Algorithms for large-scale Optimal Transport 。

### 文献阅读

文献 _A Benchmark for Discrete Optimal Transport_ 简要介绍了 OT 的提出背景和离散化过程。来源与图像相似比较或者货物运输等现实问题，一个一般的 OT 问题考察的是从 $X$ 到 $Y$ 的可测映射 $T$，希望 $T$ 具有保测度的性质，并且在给定的费用函数下使得总费用 $\int_{\Omega} C(x, T(x))d\mu$ 最小。稍简化的版本是将 $X$ 与 $Y$ 离散化，那么 $T$ 就是限制在一些离散点上的映射，因此离散化后的 OT 具有形式：

$$
min \sum_{i=1}^m \sum_{j=1}^n c_{ij} \pi_{ij} \\ s.t. \sum_{j=1}^n \pi_{ij}=\mu_i, \forall i=1, \dots m \\ \sum_{i=1}^m \pi_{ij}=\nu_j, \forall j=1, \dots n \\ \pi_{ij} \ge 0.
$$
我们可以使用解决 LP 的方法来解决离散化的 OT 问题；然而这仅仅是理论上可行的转化，由于一般图像维数太过庞大，直接使用单纯性法或者内点法是费时费空间的。另一方面的问题是离散 OT 问题有很多等式约束，但是一般化的方法没有能够很好的使用这些性质，因此我们需要使用以下的算法。

### 算法分析

- Transportation simplex，或称改进的单纯性算法，主要改进在于在迭代中选取最小的行元素，和将原问题表述为图上的传输问题来减少计算。
- Shortlist method，是一种上述方法的再改进，按照费用选定了一个较优的初始可行解。
- Shielding neighborhood method，通过解一个原问题的稀疏子集实现前几部的快速计算。
- AHA method，将2-欧氏距离的 OT 问题等价转化为一个无约束的凸优化问题。

### 课题展望

1. 将阅读参考文献的剩余部分，实现上述四种算法，并与一般解法比较性能差异。

2. 探索在给定一些的良好条件下是否能找到时间/空间复杂度较小的改进。

   ---

   目前完成了大作业要求的第一题内容，以下是代码实现与测试结果（条件 $(m,n)=(100,200)$ ）：

```matlab
clear
clc

%%
m = 100;
n = 100;

%%
% The object function f, from cost matrix c
cost = rand(m, n);
f = reshape(cost, m*n, 1);

% Equation constraint (mu - m equations; nu - n equations)
% Make sure that mu and nu have the same sum
mu = rand(m, 1);
mu = mu / sum(mu);
nu = rand(n, 1);
nu = nu / sum(nu);

% Equation coeff
Mucoeff = zeros(m, m*n);
for i = 1:m
    Mucoeff(i, i:m:m*n) = 1;
end
Nucoeff = zeros(n, m*n);
for i = 1:n
    Nucoeff(i, (i-1)*m+1:i*m) = 1;
end

%% Calling mosek with simplex method
A = -eye(m*n);
b = zeros(m*n, 1);
B = [Mucoeff; Nucoeff];
c = [mu; nu];
l = zeros(m*n, 1);
u = [];
x0 = [];
options.Simplex = 'on';
% options.Diagnostics = 'on';
options.Display = 'on';

tic
[x, fval, exitflag, output, lambda] = linprog(f, A, b, B, c, l, u, x0, options);
% Result
fprintf('Simplex method done!\n\tCost %e\n', fval);
toc

%% Calling mosek with interior point method
A = -eye(m*n);
b = zeros(m*n, 1);
B = [Mucoeff; Nucoeff];
c = [mu; nu];
l = zeros(m*n, 1);
u = [];
x0 = [];
options.Interior = 'on';
% options.Diagnostics = 'on';
options.Display = 'on';

tic
[x, fval, exitflag, output, lambda] = linprog(f, A, b, B, c, l, u, x0, options);
% Result
fprintf('Interior point method done!\n\tCost %e\n', fval);
toc
```

```plain
Simplex method done!
	Cost 2.373324e-02
Elapsed time is 0.306041 seconds.
Interior point method done!
	Cost 2.373324e-02
Elapsed time is 0.314385 seconds.
```

可见差异不是很显著。