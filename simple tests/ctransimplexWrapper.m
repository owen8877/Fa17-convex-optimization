function [x, out, val] = ctransimplexWrapper(x0, C, mu, nu, opts)
    [m, n] = size(C);
    [index, v] = ctransimplex(x0, C, mu, nu, opts);
    X = sparse(double(index(1, :))+1, double(index(2, :))+1, v);
    x = reshape(X, m*n, 1);
    out = full(sum(sum(C .* X)));
    val = [];
end