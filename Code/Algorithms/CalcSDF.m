function [sdf, contour] = CalcSDF( labeledImage, contour )
% Calculate Sign Distance Function transformation of input labeled image.
% Signature:
%   [sdf, contour] = CalcSDF( labeledImage, contour )
% Input:
%   labeledImage: labeled or BW image.
%   contour: labeled / BW contor of the input image (Optional).
% Output:
%   sdf: Sign Distance Function for the BW of the input image.
%   contour: labeled contor of the input image.
% Example:
%     I = imread('coins.png'); I = I(110:245, 55:160);
%     BW = imfill(imclearborder(im2bw(I)),'holes');
%     [phi, contour] = CalcSDF(BW);
%     f = figure; 
%     subplot(4,1,1); imagesc(I); title('Intensity Image'); colormap gray;
%     subplot(4,1,2); imagesc(BW); title('Input binary image');
%     subplot(4,1,3); imagesc(contour); title('Output contour image');
%     subplot(4,1,4); imagesc(phi); title('Output SDF'); colormap jet; colorbar
%     set(f, 'Position', get(0, 'Screensize'));
% See also: CalcContour; bwdist; im2bw; bwlabel
%
% T. Gilad, 2017
%%
    % Input contour:
    if(~exist('contour', 'var') || isempty(contour))
        contour = abs(labeledImage - imerode(labeledImage,ones(3)));
    end
    
    % Convert to binary:
    bwImage = (labeledImage ~= 0);
    bwContour = (contour ~= 0);
    
    % Distance of each pixel from contour (Euclidean):
    dist = bwdist(bwContour);
    
    % Inside contour - positive, outside - negative, on contour - zero:
    sdf = dist.*bwImage - dist.*(~bwImage);    
end

