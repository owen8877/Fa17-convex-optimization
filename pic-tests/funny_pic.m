einstein = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/einstein_resize.jpg'));
tutou = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/tutou.jpg'));
dark = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/dark.jpg'));
funny2 = double(imread('/home/xdroid/repo/convex-optimization/pic-tests/funny/funny2.jpg'));


resolution = 64;
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

mex -v CXXFLAGS='$CXXFLAGS -Wall -std=c++11 -O2' shielding.cpp
mu = reshape(dark, n, 1);
mu = mu * n / sum(mu);
nu = reshape(funny2, n, 1);
nu = nu * n / sum(nu);
[index, v] = shielding([], cost, mu, nu, []);
X = sparse(double(index(1, :))+1, double(index(2, :))+1, v);

index = double(index + 1);
for alpha = linspace(0, 1, 6)
    cor = 0.06 * min(alpha, 1-alpha)*2;
    bor = 0.12 * min(alpha, 1-alpha)*2;
    cen = 1 - 4*(cor+bor);
    conv = [-1 0 1 1 1 0 -1 -1 0; 1 1 1 0 -1 -1 -1 0 0; cor bor cor bor cor bor cor bor cen];
    figure
    inter = zeros(resolution, resolution);
    colormap gray
    for i = [double(index); v]
        [fromx, fromy] = ind2sub([resolution, resolution], i(1));
        [tox, toy] = ind2sub([resolution, resolution], i(2));
        fx = fix(fromx*(1-alpha)+tox*alpha);
        fy = fix(fromy*(1-alpha)+toy*alpha);
        for j = 1:9
            ffx = max(1, min(resolution, fx+conv(1, j)));
            ffy = max(1, min(resolution, fy+conv(2, j)));
            inter(ffx, ffy) = inter(ffx, ffy) + conv(3, j) * i(3);
        end
    end
    imagesc(inter);
end