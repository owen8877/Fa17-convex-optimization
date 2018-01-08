function [x, out, val] = transimplexWrapper(x0, C, mu, nu, opts)
    [m, n] = size(C);
    if ~isfield(opts, 'method')
        opts.method = 'normal';
    end
    if isfield(opts, 'extra')
        opts_ = opts.extra;
    end
    switch opts.method
        case 'normal'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2' normaltransimplex.cpp
            [index, v] = normaltransimplex(x0, C, mu, nu, opts_);
        case 'improved'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2' improvedtransimplex.cpp
            [index, v] = improvedtransimplex(x0, C, mu, nu, opts_);
        case 'minimalrow'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2' minimalrowtransimplex.cpp
            [index, v] = minimalrowtransimplex(x0, C, mu, nu, opts_);
        case 'shortlist'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2' shortransimplex.cpp
            [index, v] = shortransimplex(x0, C, mu, nu, opts_);
        case 'shielding'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2' shielding.cpp
            [index, v] = shielding(x0, C, mu, nu, opts_);
        case 'multiscale'
            mex CXXFLAGS='$CXXFLAGS -std=c++14 -O2' multiscale_matlab.cpp
            [index, v] = multiscale_matlab(x0, C, mu, nu, opts_);
        case 'cplex'
            mex CXXFLAGS='$CXXFLAGS -std=c++11 -O2 -m64 -fPIC -fno-strict-aliasing' ...
                -I"/opt/ibm/ILOG/CPLEX_Studio128/cplex/include" ...
                -L"/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic" ...
                -lcplex -lm -lpthread -ldl ...
                cplex_matlab.cpp
            [index, v] = cplex_matlab(x0, C, mu, nu, opts_);
    end
    X = sparse(double(index(1, :))+1, double(index(2, :))+1, v);
    x = reshape(X, m*n, 1);
    out = full(sum(sum(C .* X)));
    val = [];
end