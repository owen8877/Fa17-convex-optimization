resolutions = [16 32 64];
m = numel(resolutions);

costs = cell(m, 1);
for itr = 1:m
    resolution = resolutions(itr);
    n = resolution^2;

    cost = zeros(n, n);
    for i = 1:resolution
        for j = 1:resolution
            for k = 1:resolution
                for l = 1:resolution
                    cost(i+(j-1)*resolution, k+(l-1)*resolution) = (i-k)^2 + (j-l)^2;
                end
            end
        end
    end
    costs{itr} = cost;
end

coststruct.resolutions = resolutions;
coststruct.data = costs;

save ../data/cost.mat -struct coststruct