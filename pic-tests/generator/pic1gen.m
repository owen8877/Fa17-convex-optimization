resolutions = [16 32 64];
m = numel(resolutions);
k = 5; % pic numbers of one resolution

pic = cell(m, k);
for itr = 1:m
    resolution = resolutions(itr);
    n = resolution^2;
    for l = 1:k
        mu = rand(n, 1);
        mu = mu * n / sum(mu);
        pic{l, itr} = mu;
    end
end

picstruct.resolutions = resolutions;
picsrtuct.picN = k;
picstruct.data = pic;

save ../data/pic1.mat -struct picstruct