resolutions = [8 16 24 32 64];
m = numel(resolutions);
k = 3; % pic numbers of one resolution

pic = cell(m, k);
for itr = 1:m
    resolution = resolutions(itr);
    n = resolution^2;
    for l = 1:k
        mu = rand(n, 1);
        mu = mu / sum(mu);
        mu(end) = 1 - sum(mu(1:n-1));
        pic{l, itr} = mu;
    end
end

picstruct.resolutions = resolutions;
picstruct.picN = k;
picstruct.data = pic;

save ../data/pic1.mat -struct picstruct