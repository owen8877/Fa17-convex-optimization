function pic8_10gen(resolutions, k, fix_end, series)
    m = numel(resolutions);
    prefix = '/home/xdroid/repo/convex-optimization/pic-tests/generator/';

    pic = cell(k, m);
    for itr = 1:m
        resolution = resolutions(itr);
        n = resolution^2;
        for l = 1:k
            mu = double(imread([prefix int2str(series) '/' int2str(l) '_' int2str(resolution) '.jpg'])) + 1;
            mu = mu / sum(sum(mu));
            mu = mu + 1e-8*rand(resolution, resolution);
            mu = mu / sum(sum(mu));
            fprintf('Checking...Sum is %f.\n', sum(sum(mu)));
            pic{l, itr} = fix_end(reshape(mu, n, 1));
        end
    end

    picstruct.resolutions = resolutions;
    picstruct.picN = k;
    picstruct.data = pic;

    save(['../data/pic' int2str(series) '.mat'], '-struct', 'picstruct')
end