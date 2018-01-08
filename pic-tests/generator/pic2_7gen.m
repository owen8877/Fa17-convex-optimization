function pic2_7gen(resolutions, k, fix_end, source, dest)
    load(source)
    m = numel(resolutions);

    importFromR = cell(6, 1);
    importFromR{1} = picx8;
    importFromR{2} = picx16;
    importFromR{3} = picx24;
    importFromR{4} = picx32;
    importFromR{5} = picx64;
    importFromR{6} = picx96;
    pic = cell(k, m);
    for itr = 1:m
        resolution = resolutions(itr);
        tmp = importFromR{itr};
        n = resolution^2;
        for l = 1:k
            mu = tmp((l-1)*resolution+1:l*resolution, :);
            mu = mu + 1e-10;
            mu = mu / sum(sum(mu));
            fprintf('Checking...Sum is %f.\n', sum(sum(mu)));
            pic{l, itr} = fix_end(reshape(mu, n, 1));
        end
    end

    picstruct.resolutions = resolutions;
    picstruct.picN = k;
    picstruct.data = pic;

    save(dest, '-struct', 'picstruct')
end