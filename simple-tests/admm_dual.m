function [x, out, val] = admm_dual(x0, cost, mu, nu, opts)
    opts = default(opts);
    
    [m, n] = size(cost);
    mpn = m+n; mtn = m*n;
    f = reshape(cost, m*n, 1);
    mui = repmat(1:m, 1, n);
    muj = 1:m*n;
    muval = ones(1, m*n);
    muCoeff = sparse(mui, muj, muval);
    nui = reshape(repmat(1:n, m, 1), 1, m*n);
    nuj = 1:m*n;
    nuval = ones(1, m*n);
    nuCoeff = sparse(nui, nuj, nuval);
    T = [muCoeff; nuCoeff];
    s = [mu; nu];
    
    x = x0; hx = x;
    itr = 1;
    z = zeros(mtn, 1); hz = z;
    y = zeros(mpn, 1); hy = y;
    beta = 0.05;
    alpha = 1;
    eps = 1e-9;
    val = [];
    
    while true % && itr < 10000
        if ~opts.nesterov
            lx = x;
            ly = y;
            lz = z;
            for j = 1:1
                % update y
                b = (T*x-s)/beta-T*(f+lz);
%                 y = cg(ly, @(v) T*(T'*v) , b);
                y = cg(ly, @(v) T*(T'*v) + eps*v, b);
%                 sy = (sum(b(m+1:end))*m-sum(b(1:m))) / (m*n-1);
%                 sz = (sum(b(1:m)*n)-sum(b(m+1:m+n))) / (m*n-1);
%                 y = [b(1:m)-sz; b(m+1:m+n)-sy];
                % update z
                z = min(x/beta-(T'*y+f), 0);
            end
            % update x
            x = x - beta*(z+f+T'*y);
        else
            lx = x;
            ly = y;
            lz = z;
            % update y
            b = (T*hx-s)*beta-T*(f+hz);
            y = cg(ly, @(v) T*(T'*v) + eps*v, b);
            % update z
            z = min(-hx, 0);
            % update x
            x = hx + (z-f-T'*y)/beta;
            alpha_ = (1+sqrt(1+4*alpha^2))/2;
            hz = z + (alpha-1)/alpha_ * (z-lz);
            hx = x + (alpha-1)/alpha_ * (x-lx);
            alpha = alpha_;
%             lx = x;
%             ly = y;
%             lz = z;
%             for j = 1:1
%                 % update z
%                 z = min(hx*beta-(T'*hy+f), 0);
%                 % update y
%                 b = (T*hx-s)*beta-T*(f+z);
%                 y = cg(ly, @(v) T*(T'*v) + eps*v, b);
%             end
%             % update x
%             x = hx - (z+f+T'*y)/beta; 
%             alpha_ = (1+sqrt(1+4*alpha^2))/2;
%             hy = y + (alpha-1)/alpha_ * (y-ly);
%             hx = x + (alpha-1)/alpha_ * (x-lx);
%             alpha = alpha_;
        end
        
        if mod(itr, 100) == 0
            out = f'*x;
            constraintErr = norm(s-T*x);
            fprintf('%d\t%.8e\t%.1e\t%.2f\t%.2e\n', ...
                itr, out, constraintErr, sum(x>=0)/mtn, alpha);
            if constraintErr < opts.tor
                break
            end
        end
        
        itr = itr + 1;
    end
end

function opts = default(opts)
    if ~isfield(opts, 'nesterov')
        opts.nesterov = false;
    end
    if ~isfield(opts, 'tor')
        opts.tor = 1e-5;
    end
end