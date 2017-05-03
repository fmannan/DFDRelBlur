% global_affine_flow_S1T2.m
% Fahim Mannan (fmannan@cim.mcgill.ca)
function T = global_affine_flow_S1T2(im1, im2, params)
% 1 Scale (equal scaling in both dimension) and 2 Translation parameters
% same as global_affine_flow except entries corresponding to a2 and a4 are
% removed
im1 = mean(im1, 3);
im2 = mean(im2, 3);

[Ix, Iy] = gradient(im1);
It = im2 - im1;

Ix2 = Ix .* Ix;
Iy2 = Iy .* Iy;
IxIy = Ix .* Iy;

[X, Y] = meshgrid(1:size(im1, 2), 1:size(im1, 1));

sIx2 = nansum(Ix2(:));
sIx2_x = nansum(Ix2(:) .* X(:));
sIx2_x2 = nansum(Ix2(:) .* X(:).^2);
%sIx2_y = nansum(Ix2(:) .* Y(:));
%sIx2_y2 = nansum(Ix2(:) .* Y(:).^2);
%sIx2_xy = nansum(Ix2(:) .* X(:) .* Y(:));

sIy2 = nansum(Iy2(:));
%sIy2_x = nansum(Iy2(:) .* X(:));
%sIy2_x2 = nansum(Iy2(:) .* X(:).^2);
sIy2_y = nansum(Iy2(:) .* Y(:));
sIy2_y2 = nansum(Iy2(:) .* Y(:).^2);
%sIy2_xy = nansum(Iy2(:) .* X(:) .* Y(:));

sIxIy = nansum(IxIy(:));
sIxIy_x = nansum(IxIy(:) .* X(:));
%sIxIy_x2 = nansum(IxIy(:) .* X(:).^2);
sIxIy_y = nansum(IxIy(:) .* Y(:));
%sIxIy_y2 = nansum(IxIy(:) .* Y(:).^2);
sIxIy_xy = nansum(IxIy(:) .* X(:) .* Y(:));

A = [sIx2 , sIx2_x + sIxIy_y,  sIxIy    ;
    sIx2_x + sIxIy_y , sIx2_x2 + 2 * sIxIy_xy + sIy2_y2,  sIxIy_x + sIy2_y ;
    sIxIy , sIxIy_x + sIy2_y ,  sIy2  ];
    
b = -[nansum(Ix(:) .* It(:)) ; 
      nansum(Ix(:) .* It(:) .* X(:)) + nansum(Iy(:) .* It(:) .* Y(:)); 
      nansum(Iy(:) .* It(:)) ];
  
a = A \ b;

T = [(1 + a(2)) 0 a(1); 0 (1 + a(2)) a(3); 0 0 1];

