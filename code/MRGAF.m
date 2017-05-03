% MRGAF.m
% Fahim Mannan (fmannan@cim.mcgill.ca)
function [Afinal, im1final, im2final] = MRGAF(im1_orig, im2_orig, params)
% Multiresolution Global Affine Flow
nLevels = params.nLevels;
scale = params.scale;
innerloop_threshold = 1e-6;
maxIter = 10;
if(isfield(params, 'InnerLoopThreshold'))
    innerloop_threshold = params.InnerLoopThreshold;
end

if(isfield(params, 'MaxInnerIteration'))
    maxIter = params.MaxInnerIteration;
end

algName = 'GLOBAL_AFFINE';
if(isfield(params, 'ALGORITHM'))
    algName = params.ALGORITHM
end
if(isfield(params, 'AInit') && ~isempty(params.AInit))
    A0 = params.AInit;
    [X, Y] = meshgrid(1:size(im1_orig, 2), 1:size(im1_orig, 1));
    Pt_old = [X(:) Y(:)]';
    Pt_old(3,:) = 1;
    Pt_new = A0 \ Pt_old;

    Xnew = reshape(Pt_new(1,:) ./ Pt_new(3,:), [size(im1_orig, 1), size(im1_orig, 2)]);
    Ynew = reshape(Pt_new(2,:) ./ Pt_new(3,:), [size(im1_orig, 1), size(im1_orig, 2)]);
    im1warp = nan(size(im1_orig));
    for k = 1:3
        im1warp(:, :, k) = interp2(X, Y, double(im1_orig(:, :, k)), Xnew, Ynew);
    end
    border = sum(im1warp, 3);
    [h, w] = size(border);
    startWidth = max(sum(isnan(border(ceil(h/2), 1:ceil(w/2))), 2));
    endWidth = max(sum(isnan(border(ceil(h/2), ceil(w/2):end)), 2));
    startHeight = max(sum(isnan(border(1:ceil(h/2), ceil(w/2))), 1));
    endHeight = max(sum(isnan(border(ceil(h/2):end, ceil(w/2))), 1));
        
    im1_orig = double(im1warp(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
    im2_orig = double(im2_orig(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
else
    A0 = eye(3, 3);
end
Afinal = eye(3, 3);

pyrIm1 = cell(1, nLevels + 1);
pyrIm2 = cell(1, nLevels + 1);

pyrIm1{1} = im1_orig;
pyrIm2{1} = im2_orig;
for l = 1:nLevels
    pyrIm1{l + 1} = imresize(pyrIm1{l}, scale);
    pyrIm2{l + 1} = imresize(pyrIm2{l}, scale);
end

for l = nLevels+1:-1:1
    im1 = pyrIm1{l};
    im2 = pyrIm2{l};

    if(l <= nLevels && nLevels > 0)
        % upscale translation
        Afinal(1:2, 3) = Afinal(1:2, 3) * 1/scale;
        [X, Y] = meshgrid(1:size(im1, 2), 1:size(im1, 1));
        Pt_old = [X(:) Y(:)]';
        Pt_old(3,:) = 1;
        Pt_new = Afinal \ Pt_old;

        Xnew = reshape(Pt_new(1,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
        Ynew = reshape(Pt_new(2,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
        im1warp = nan(size(im1));
        for k = 1:3
            im1warp(:, :, k) = interp2(X, Y, double(im1(:, :, k)), Xnew, Ynew);
        end
        %im1 = im1warp;
        border = sum(im1warp, 3);
        [h, w] = size(border);
        startWidth = max(sum(isnan(border(ceil(h/2), 1:ceil(w/2))), 2));
        endWidth = max(sum(isnan(border(ceil(h/2), ceil(w/2):end)), 2));
        startHeight = max(sum(isnan(border(1:ceil(h/2), ceil(w/2))), 1));
        endHeight = max(sum(isnan(border(ceil(h/2):end, ceil(w/2))), 1));

        im1 = double(im1warp(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
        im2 = double(im2(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
    end
    for i = 1:maxIter
        if(strcmpi(algName, 'GLOBAL_AFFINE_S2T2'))
            A = global_affine_flow_S2T2(im1, im2);
        elseif(strcmpi(algName, 'GLOBAL_AFFINE_S1T2'))
            A = global_affine_flow_S1T2(im1, im2);
        elseif(strcmpi(algName, 'GLOBAL_AFFINE_S1'))
            A = global_affine_flow_S1(im1, im2);
        else
            A = global_affine_flow(im1, im2); % pyrIm2{l});
        end

        A(3,:) = [ 0, 0 , 1];
        %A
        %norm(reshape(A - eye(3, 3), 1, []), 2)
        err = norm(reshape(A - eye(3, 3), 1, []), 1);
        if(err < innerloop_threshold)
            break;
        end
        Afinal = A * Afinal; % the final transformation for the current level

        [X, Y] = meshgrid(1:size(im1, 2), 1:size(im1, 1));
        Pt_old = [X(:) Y(:)]';
        Pt_old(3,:) = 1;
        Pt_new = A \ Pt_old;

        Xnew = reshape(Pt_new(1,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
        Ynew = reshape(Pt_new(2,:) ./ Pt_new(3,:), [size(im1, 1), size(im1, 2)]);
        im1warp = nan(size(im1));
        for k = 1:3
            im1warp(:, :, k) = interp2(X, Y, double(im1(:, :, k)), Xnew, Ynew);
        end
        %figure
        %imshow(uint8(im1warp));

        %im1 = im1warp;
        border = sum(im1warp, 3);
        [h, w] = size(border);
        startWidth = max(sum(isnan(border(ceil(h/2), 1:ceil(w/2))), 2));
        endWidth = max(sum(isnan(border(ceil(h/2), ceil(w/2):end)), 2));
        startHeight = max(sum(isnan(border(1:ceil(h/2), ceil(w/2))), 1));
        endHeight = max(sum(isnan(border(ceil(h/2):end, ceil(w/2))), 1));

        im1 = double(im1warp(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
        im2 = double(im2(startHeight+1:end-endHeight-1, startWidth + 1:end-endWidth-1, :)); %im1warp;
    end
    im1final = im1;
    im2final = im2;
end
Afinal = Afinal * A0;
