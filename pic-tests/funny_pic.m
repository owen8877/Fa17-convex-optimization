% einstein = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/einstein_resize.jpg'));
% tutou = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/tutou.jpg'));
% dark = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/dark.jpg'));
% dawei = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/funny2.jpg'));
% van = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/van.jpg'));
% bill = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/bill.jpg'));
% cup1 = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/generator/10/4_64.jpg'));
% cup2 = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/generator/10/5_64.jpg'));
frame1 = double(imread('/home/xdroid/Downloads/a/img1.png'));
frame2 = double(imread('/home/xdroid/Downloads/a/img2.png'));

% resolution = 64;
% n = resolution^2;
% cost = zeros(n, n);
% for i = 1:resolution
%     for j = 1:resolution
%         for k = 1:resolution
%             for l = 1:resolution
%                 cost(i+(j-1)*resolution, k+(l-1)*resolution) = (i-k)^2 + (j-l)^2;
%             end
%         end
%     end
% end
width = 113; height = 117;
m = width * height; n = m;
cost = zeros(m, n);
for i = 1:width
    for j = 1:width
        for k = 1:height
            for l = 1:height
                cost((i-1)*height+k, (j-1)*height+l) = (i-j)^2 + (k-l)^2;
            end
        end
    end
end
mex CXXFLAGS='$CXXFLAGS -std=c++14 -Ofast' shielding.cpp
mu = reshape(frame1, n, 1);
mu = mu / sum(mu);
nu = reshape(frame2, n, 1);
nu = nu / sum(nu);
[index, v] = shielding([], cost, mu, nu, false, height);
X = sparse(double(index(1, :))+1, double(index(2, :))+1, v);
index = double(index + 1);

% all = [];
% for alpha = linspace(0, 1, 6)
%     cor = 0.00 * min(alpha, 1-alpha)*2;
%     bor = 0.00 * min(alpha, 1-alpha)*2;
%     cen = 1 - 4*(cor+bor);
%     conv = [-1 0 1 1 1 0 -1 -1 0; 1 1 1 0 -1 -1 -1 0 0; cor bor cor bor cor bor cor bor cen];
%     inter = zeros(resolution, resolution);
%     for i = [double(index); v]
%         [fromx, fromy] = ind2sub([resolution, resolution], i(1));
%         [tox, toy] = ind2sub([resolution, resolution], i(2));
%         fx = fix(fromx*(1-alpha)+tox*alpha);
%         fy = fix(fromy*(1-alpha)+toy*alpha);
%         for j = 1:9
%             ffx = max(1, min(resolution, fx+conv(1, j)));
%             ffy = max(1, min(resolution, fy+conv(2, j)));
%             inter(ffx, ffy) = inter(ffx, ffy) + conv(3, j) * i(3);
%         end
%     end
% %     imagesc(inter);
%     all = [all inter];
% end
% figure
% colormap gray
% imagesc(all);
% axis image
% brighten(0.4)