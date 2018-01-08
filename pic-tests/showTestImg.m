resolution = 96;
allImage = zeros(96*5, 96*10);

correction = [1, 1, 1, 1, 1, 1, 1, 1.5, 3, 3];
for testCase = 1:10
    picData = load(['data/pic', int2str(testCase), '.mat']);
    for index = 1:5
        allImage((index-1)*resolution+1:index*resolution, (testCase-1)*resolution+1:testCase*resolution) = reshape(picData.data{index, 6}, resolution, resolution) * correction(testCase);
    end
end

colormap gray
imagesc(allImage);
brighten(0.1)
axis image