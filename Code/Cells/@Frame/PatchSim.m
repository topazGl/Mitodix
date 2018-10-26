function sim = PatchSim( im1, mask1, im2, mask2, registerFlag, r )
% Calculate weighted similarity between 2 patches.
%   sim = PatchSim( im1, mask1, im2, mask2 )
% Inputs:
%   im1, im2 - patches (gray scale).
%   mask1, mask2 - binary mask images ('1' = cell pixel).
%   registerFlag = If '1': register images first according to their fuzzy masked
%   version with rigid transformation. '2': registration using translation. '0': No registration (Default=0).
% Outputs:
%   weighted mean ratio correlation similarity.
    if(~exist('registerFlag', 'var') || isempty(registerFlag))
        registerFlag = 0;
    end
    if(~exist('r', 'var') || isempty(r))
        r = 1;
    end
    
    if(registerFlag)
        if(registerFlag == 2)
            [im2, mask2] = RegisterPatches(im1, im2, mask1, mask2, 'translation');
        else
            [im2, mask2] = RegisterPatches(im1, im2, mask1, mask2); 
        end
    end
    mask = (mask1 | mask2);
    [~,~,~,weights] = Frame.TransformLabeledImage(mask, 1);
    [wcorr,~,~,sim] = wcorr2( im1, im2, weights, r, 0 );
%     if(abs(wcorr-sim) > 0.2)
%         [wcorr1,~,~,sim1] = wcorr2( im1, im2, weights, 1 );
%         corr = corr2( im1, im2 );
%         wcorrm = wcorr2( im1, im2, weights,[], 0);
%         M = max([im1(:); im2(:)]);
%         disp(['corr=', num2str(corr),'; wcorrm=', num2str(wcorrm), '; wcorr0.3=', num2str(wcorr), '; sim0.3=', num2str(sim),'; wcorr1=', num2str(wcorr1), '; sim1=', num2str(sim1)]);
%         f1 = figure; subplot(2,1,1); imagesc(weights/max(weights(:))); colorbar; ax1 = gca; subplot(2,1,2); imshowpair(im1.*weights, im2.*weights); ax2 = gca; linkaxes([ax1, ax2]);
%         f2 = figure; subplot(2,1,1); imagesc(im1/M); cc = get(gca, 'clim'); ax3 = gca; colorbar; subplot(2,1,2); imagesc(im2/M); set(gca, 'clim', cc); colorbar; linkaxes([ax3, gca]); colormap gray;
%         close(f1); close(f2);
%     end
end

function [im, mask] = RegisterPatches(imFixed, imMoving, maskFixed, maskMoving, transType)
    if(~exist('transType', 'var') || isempty(transType))
        transType = 'rigid';
    end
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = 1000;
    optimizer.MaximumStepLength = 10^-3;
    se = strel('disk',1);
    
    [~,~,~,weightsFixed] = Frame.TransformLabeledImage(maskFixed, 1);
    [~,~,~,weightsMoving] = Frame.TransformLabeledImage(maskMoving, 1);
    maskedFixed = double(imFixed).*weightsFixed/max(weightsFixed(:));
    maskedMoving = double(imMoving).*weightsMoving/max(weightsMoving(:));

    Rmoving = imref2d(size(maskedFixed));
    Rfixed = imref2d(size(maskedMoving));
    tform  = imregtform(maskedMoving, maskedFixed, transType, optimizer, metric);

    im = imwarp(imMoving, Rmoving, tform, 'cubic', 'OutputView',Rfixed );
    im(im < 0) = 0;
    mask = imwarp(maskMoving, Rmoving, tform, 'cubic', 'OutputView',Rfixed );
    mask = imclose((mask ~= 0), se);
end
