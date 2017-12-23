function [x, out, val] = auglag(x0, cost, mu, nu, opts)
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
    
    x = x0; lx = x; hx = x;
    itr = 1;
    z = max(x, 0); lz = z; hz = z;
    y = zeros(mpn+mtn, 1); ly = y; hy = y;
    alpha = 1;
    beta = 1e-3;
    eta = 0.999;
    oldcrit = Inf;
    restartTimes = 0;
    val = zeros(3, 15000);
    
    while true && itr < 1000000
        if ~opts.nesterov
            for j = 1:5
                % update x
                v1 = s + y(1:mpn)/beta;
                v2 = z + y(mpn+1:end)/beta;
                nx = cg(x, @(vt) T'*(T*vt) + vt, T'*v1 + v2 - f/beta);
                % update z
                z = max(0, x-y(mpn+1:end)/beta);
                x = nx;
            end
            % update y
            y = y + [s-T*x; z-x];
        else
%             % update x
%             v1 = s + hy(1:mpn)/beta;
%             v2 = hz + hy(mpn+1:mpn+mtn)/beta;
%             x = cg(x, @(vt) T'*(T*vt) + vt, T'*v1 + v2 - f/beta);
%             % update z
%             z = max(0, x-hy(mpn+1:mpn+mtn)/beta);
%             % update y
%             y = hy + beta*[s-T*x; z-x];
%             crit = (norm(y-hy)^2)/beta + beta*norm(z-hz)^2;
%             if crit < eta*oldcrit
%                 alpha_ = (1+sqrt(1+4*alpha^2))/2;
%                 hz = z + (alpha-1)/alpha_*(z-lz);
%                 hy = y + (alpha-1)/alpha_*(y-ly);
%                 lz = z;
%                 ly = y;
%                 alpha = alpha_;
%                 oldcrit = crit;
%             else
%                 alpha = 1;
%                 hz = lz;
%                 lz = z;
%                 hy = ly;
%                 ly = y;
%                 oldcrit = oldcrit / eta;
%                 % fprintf('Restarted!\n')
%                 if restartTimes < 100
%                     restartTimes = restartTimes + 1;
%                 else
%                     break
%                 end
%             end
            % update z
            z = max(0, hx-hy(mpn+1:mpn+mtn)/beta);
            % update x
            v1 = s + hy(1:mpn)/beta;
            v2 = z + hy(mpn+1:mpn+mtn)/beta;
            x = cg(lx, @(vt) T'*(T*vt) + vt, T'*v1 + v2 - f/beta);
            % x = pcg(@(vt) T'*(T*vt) + vt, T'*v1 + v2 - f/beta, 1e-10, [], [], [], x);
            % update y
            y = hy + 1*[s-T*x; z-x];
            crit = norm(y-hy)^2/beta + beta*norm(x-hx)^2;
            if crit < eta*oldcrit
                alpha_ = (1+sqrt(1+4*alpha^2))/2;
                hx = x + (alpha-1)/alpha_*(x-lx);
                hy = y + (alpha-1)/alpha_*(y-ly);
                lx = x;
                ly = y;
                alpha = alpha_;
                oldcrit = crit;
            else
                alpha = 1;
                hx = lx;
                lx = x;
                hy = ly;
                ly = y;
                oldcrit = oldcrit / eta;
                fprintf('R')
                if restartTimes < 20 || true
                    restartTimes = restartTimes + 1;
                else
                    break
                end
            end
        end
        
        if mod(itr, 100) == 0
            out = f' * x;
            fprintf('%d\t%.8e\t%.1e\t%.2f\t\n', ...
                itr, out, norm([T*x-s; x-z]), sum(x>1e-9)/mtn);
            if norm([T*x-s; x-z]) < opts.tor
                break
            end
        end
        
        % val(:, itr) = [min(1, f'*x); norm(opts.mosek-x); norm(T*x-s)];
        itr = itr + 1;
    end
    
    val = val(:, 1:itr);
    out = f'*x;
end

function opts = default(opts)
    if ~isfield(opts, 'nesterov')
        opts.nesterov = true;
    end
end