clc; clear;
resolutions = [8 16 24 32 64 96];
k = 5;

% pic1gen(resolutions, k, @fix_end);
% pic2_7gen(resolutions, k, @fix_end, 'GRF_2.mat', '../data/pic2.mat');
% pic2_7gen(resolutions, k, @fix_end, 'GRF_3.mat', '../data/pic3.mat');
% pic2_7gen(resolutions, k, @fix_end, 'GRF_4.mat', '../data/pic4.mat');
% pic2_7gen(resolutions, k, @fix_end, 'GRF_5.mat', '../data/pic5.mat');
% pic2_7gen(resolutions, k, @fix_end, 'GRF_6.mat', '../data/pic6.mat');
% pic2_7gen(resolutions, k, @fix_end, 'GRF_7.mat', '../data/pic7.mat');
pic8_10gen(resolutions, k, @fix_end, 8);
pic8_10gen(resolutions, k, @fix_end, 9);
pic8_10gen(resolutions, k, @fix_end, 10);

function v = fix_end(v0)
    v = v0;
    v(end) = 1-sum(v(1:numel(v)-1));
end