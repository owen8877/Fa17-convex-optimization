function [x, out, val] = ctransimplexWrapper(x0, C, mu, nu, opts)
    [m, n] = size(C);
    if ~isfield(opts, 'method')
        opts.method = 'normal';
    end
    switch opts.method
        case 'normal'
            mex -v CXXFLAGS='$CXXFLAGS -Wall -std=c++11 -O2' ctransimplex.cpp
            [index, v] = ctransimplex(x0, C, mu, nu, opts);
        case 'improved'
            mex -v CXXFLAGS='$CXXFLAGS -Wall -std=c++11 -O2' improvedtransimplex.cpp
            [index, v] = improvedtransimplex(x0, C, mu, nu, opts);
        case 'shortlist'
            mex -v CXXFLAGS='$CXXFLAGS -Wall -std=c++11 -O2' shortransimplex.cpp
            [index, v] = shortransimplex(x0, C, mu, nu, opts);
    end
    X = sparse(double(index(1, :))+1, double(index(2, :))+1, v);
    x = reshape(X, m*n, 1);
    out = full(sum(sum(C .* X)));
    val = [];
end