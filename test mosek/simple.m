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