% imwarp.m
% Fahim Mannan (fmannan@cim.mcgill.ca)
function im1warp = imwarp(im1, A)

[X, Y] = meshgrid(1:size(im1, 2), 1:size(im1, 1));
Pt_old = [X(:) Y(:)]';
Pt_old(3,:) = 1;
Pt_new = A \ Pt_old;

Xnew = reshape(Pt_new(1,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
Ynew = reshape(Pt_new(2,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
im1warp = nan(size(im1));
K = size(im1, 3);
for k = 1:K
    im1warp(:, :, k) = interp2(X, Y, double(im1(:, :, k)), Xnew, Ynew);
end
