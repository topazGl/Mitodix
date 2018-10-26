function [imContour, imSdf, imPhi, imSoftMask, biggerMask] = TransformLabeledImage(labeled, epsilon, mask, softThr, mskThr)
% Calculate contour, SDF, phi and soft mask for the input labeled image.
%   [imContour, imSdf, imPhi, imSoftMask, biggerMask] = TransformLabeledImage(labeled, epsilon, mask, softThr, mskThr)
% Inputs:
%   labeled = labeled image.
%   epsilon = (optional). Factor for sigmoid. 
%   mask = (optional) mask image.
%   softThr = (Default = 0). Soft mask threshold.
%   mskThr = (Default = 0). Soft mask to mask threshold.
% Ouputs:
%   imContour = image contour.
%   imSdf = image sign-distane-function transformation.
%   imPhi = image sigmoid transformation.
%   imSoftMask = image soft-mask.
%   biggerMask = binary image according to threshold on imSoftMask.
    [imSdf, imContour] = CalcSDF(labeled);
        if((exist('epsilon', 'var')) && (~isempty(epsilon)))
            imPhi = CalcPHI(imSdf, epsilon);
        else
            imPhi = CalcPHI(imSdf);
        end

        % --- calculate soft mask: ---
        imSoftMask = imPhi;
        if(~exist('softThr', 'var') || isempty(softThr))
            softThr = 0;
        end

        if(~exist('mskThr', 'var') || isempty(mskThr))
            mskThr = softThr;
        end
        biggerMask = (imSoftMask >= mskThr);

        imSoftMask(imSoftMask < softThr) = 0;

        if((nargin > 2) && (~isempty(mask)))
            % mask phi:
            imSoftMask(~mask) = 0;
        end

        imSoftMask = imSoftMask/(sum(imSoftMask(:)));
end

