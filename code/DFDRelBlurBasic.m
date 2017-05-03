% DFDRelBlurBasic.m
% For personal, academic, or educational use only
% Fahim Mannan (fmannan@cim.mcgill.ca)
% Usage:
% This script loads an image, aligns it if bAlignment = 1, followed by 
% appropriately cropping the border. Then it smooths the images if 
% bNoisePrefiltering = 1. DFD is performed by blurring both images
% isotropicially and then taking pixel-wise differences. This results in a
% cost matrix. If bCostSmoothing = 1 then this matrix is smoothed to remove
% noisy cost (this is similar to using an aggregation window). Finally the
% resulting relative blur assignment is further smoothed using the 
% PostSmoothingKernel.
% 
close all
clear
clc

%addpath '../OpticalFlow'

% Processing parameters
bAlignment = 0;             % set to 1 for image alignment 
bNoisePrefiltering = 0;
bCostSmoothing = 1;
PostSmoothingKernel = fspecial('gaussian', 7, 1);
ScaleImage = 1; %
etime = [];

% Input file specification
DirRoot = '../data/';
filenamePrefix = 'breakfast'; %'img1'; %'img4';  %'img3';
filenameSuffix = '.ppm'; %'.JPG'; %
NearSuffix = '_near'; % '_1' %
FarSuffix = '_far';   % '_inf' %
imNear = im2double(imread([DirRoot filenamePrefix NearSuffix filenameSuffix]));
imFar = im2double(imread([DirRoot filenamePrefix FarSuffix filenameSuffix]));

imNear = imresize(imNear, ScaleImage);
imFar = imresize(imFar, ScaleImage);

params.nLevels = 3;
params.scale = .5;
params.InnerLoopThreshold = 5e-4;
params.ALGORITHM = 'GLOBAL_AFFINE_S1T2';

im11 = imNear;
if(bAlignment)
    tic
    A = MRGAF(imNear, imFar, params)
    etime(1) = toc;
    display(['Global Affine Flow time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);

    tic
    im11 = imwarp(imNear, A);
    etime(2) = toc;
    display(['Image Warping time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);
end
figure
imagesc(im11);
figure
imagesc(imFar);

% Crop border
border = sum(im11, 3);
[R, C] = size(border);
startWidth = max(sum(isnan(border(ceil(R/2), 1:ceil(C/2))), 2));
endWidth = max(sum(isnan(border(ceil(R/2), ceil(C/2):end)), 2));
startHeight = max(sum(isnan(border(1:ceil(R/2), ceil(C/2))), 1));
endHeight = max(sum(isnan(border(ceil(R/2):end, ceil(C/2))), 1));

imNearC = im11(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :);
imFarC = imFar(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :);

figure
imagesc(imNearC);
figure
imagesc(imFarC);

%% Use isotropic heat equation
%result = DFD2F2DBCC(params)
I1 = imNearC;
I2 = imFarC;
dT = 0.05;
nIter = 16;

if(bNoisePrefiltering)
    tic
    NoiseFilter = fspecial('gaussian', 3, 0.5);
    I1 = imfilter(I1, NoiseFilter, 'same');
    I2 = imfilter(I2, NoiseFilter, 'same');
    etime(end + 1) = toc;
    display(['Noise Prefiltering time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);
end

tic
[imb1, imbs1] = isotropic_diffusion(I1(:,:,2), 1, dT, nIter);
etime(end + 1) = toc;
display(['Isotropic Diffusion (near image) time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);
tic
[imb2, imbs2] = isotropic_diffusion(I2(:,:,2), 1, dT, nIter);
etime(end + 1) = toc;
display(['Isotropic Diffusion (far image) time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);

%%
nSZ = size(imbs1, 3);
sigval = sqrt(2* (0:nSZ-1) * dT); % relative blur

sigVal = [-sigval(end:-1:2) sigval];

tic
errF = abs(bsxfun(@minus, imbs1(:,:,end:-1:1), double(I2(:,:,1))));
errF(:,:,nSZ:2*nSZ-1) = abs(bsxfun(@minus, imbs2, double(I1(:,:,1))));
etime(end + 1) = toc;
display(['Cost Function time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);

errFs = errF;
if(bCostSmoothing)
    H = fspecial('disk', 3);
    tic
    for i = 1:size(errF, 3)
        errFs(:,:,i) = conv2(errF(:,:,i), H, 'same');
    end
    etime(end + 1) = toc;
    display(['Cost Smoothing time = ' num2str(etime(end)) ' , Total = ' num2str(sum(etime))]);
end
[val idx] = min(errF, [], 3);
[val idx1] = min(errFs, [], 3);

RelBlur = sigVal(idx);
figure
imagesc(RelBlur)
colormap gray
title('Without Smoothing');

if(bCostSmoothing)
    RelBlurCS = sigVal(idx1);
    figure
    imagesc(RelBlurCS)
    colormap gray
    title('With Basic Smoothing');
end

if(~isempty(PostSmoothingKernel))
    RelBlurCS_PS = imfilter(RelBlurCS, PostSmoothingKernel);
    figure
    imagesc(RelBlurCS_PS)
    colormap gray
    title('With Basic Smoothing (post smoothing)');
end