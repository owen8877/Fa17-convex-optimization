function [x, out, val] = dir_mosek(~, C, mu, nu, opts)
    [m, n] = size(C);
    f = reshape(C, m*n, 1);
    
    mui = repmat(1:m, 1, n);
    muj = 1:m*n;
    muval = ones(1, m*n);
    muCoeff = sparse(mui, muj, muval);
    nui = reshape(repmat(1:n, m, 1), 1, m*n);
    nuj = 1:m*n;
    nuval = ones(1, m*n);
    nuCoeff = sparse(nui, nuj, nuval);
    
    A = [];
    b = [];
    B = [muCoeff; nuCoeff];
    c = [mu; nu];
    l = zeros(m*n, 1);
    u = [];
    x0 = [];
    switch opts.method
        case 'interior'
            options.Interior = 'on';
        case 'simplex'
            options.Simplex = 'on';
        otherwise
            options.Interior = 'on';
    end
    options.Diagnostics = 'on';
    options.Display = 'iter';

    [x, ~, ~, ~, ~] = linprog(f, A, b, B, c, l, u, x0, options);

    out = sum(sum(C.*reshape(x, m, n)));
    val = [];
end

