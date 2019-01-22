load example.mat

% all = [];
% for alpha = linspace(0, 1, 6)
%     cor = 0.05 * min(alpha, 1-alpha)*2;
%     bor = 0.10 * min(alpha, 1-alpha)*2;
%     cen = 1 - 4*(cor+bor);
%     conv = [-1 0 1 1 1 0 -1 -1 0; 1 1 1 0 -1 -1 -1 0 0; cor bor cor bor cor bor cor bor cen];
%     inter = zeros(height, width);
%     for i = [double(index); v]
%         [fromx, fromy] = ind2sub([height, width], i(1));
%         [tox, toy] = ind2sub([height, width], i(2));
%         fx = fix(fromx*(1-alpha)+tox*alpha);
%         fy = fix(fromy*(1-alpha)+toy*alpha);
%         for j = 1:9
%             ffx = max(1, min(height, fx+conv(1, j)));
%             ffy = max(1, min(width, fy+conv(2, j)));
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

[fromx, fromy] = ind2sub([height, width], index(1, :));
[tox, toy] = ind2sub([height, width], index(2, :));
delta = normC([tox - fromx; toy - fromy]);
quiver(fromx, fromy, delta(1, :).*v, delta(2, :).*v);


function AA = normC(A)
    AA = A./repmat(sqrt(sum(A.^2,1)),size(A,1),1);
end