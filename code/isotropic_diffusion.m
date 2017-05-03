% isotropic_diffusion.m
% Fahim Mannan (fmannan@cim.mcgill.ca)
function [blurOut, blurredIm]= isotropic_diffusion(im, c, dT, N)
% usage : blurIm = isotropic_diffusion(imd, 1, .1, 5);
blurredIm = zeros(size(im, 1) + 2, size(im, 2) + 2, N + 1);
blurredIm(:,:,1) = padarray(im, [1, 1], 'replicate'); 
%h = figure;
for iter = 1:N
    [uxx, uyy] = central_diff2(blurredIm(:,:,iter));
    tmp = blurredIm(:,:,iter) + ...
        dT * c * (uxx + uyy);
    blurredIm(:,:,iter+1) = padarray(tmp(2:end-1,2:end-1,:), [1, 1], 'replicate'); % no-flux boundary condition
%     figure(h)
%     imagesc(blurredIm(2:end-1,2:end-1,iter+1))
%     axis image
%     colormap gray
%     title(['Iteration ' num2str(iter) ', T =' num2str(iter * dT) ', \sigma=' num2str(sqrt(2*iter*dT))]);
end
blurredIm = blurredIm(2:end-1,2:end-1,:);
blurOut = blurredIm(:,:,N+1);
