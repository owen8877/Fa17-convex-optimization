function pic1gen(resolutions, k, fix_end)
    m = numel(resolutions);

    pic = cell(k, m);
    for itr = 1:m
        resolution = resolutions(itr);
        n = resolution^2;
        for l = 1:k
            mu = rand(n, 1);
            mu = mu / sum(mu);
            mu = fix_end(mu);
            pic{l, itr} = mu;
        end
    end

    picstruct.resolutions = resolutions;
    picstruct.picN = k;
    picstruct.data = pic;

    save ../data/pic1.mat -struct picstruct
end