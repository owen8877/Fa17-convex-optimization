function x = cg(x0, A, b)
    itr = 1;
    x = x0;
    r = b - A(x);
    lr = r;
    while norm(r) > 1e-20
        if itr == 1
            p = r;
        else
            p = r + (r'*r)/(lr'*lr)*p;
        end
        alpha = (r'*r) / (p'*A(p));
        x = x + alpha*p;
        lr = r;
        r = lr - alpha*A(p);
        itr = itr + 1;
    end
end

