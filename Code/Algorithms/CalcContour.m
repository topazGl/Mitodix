function contour = CalcContour( labeledImage )
% Calculate conour of input labeled image.
% Signature:
%   contour = CalcContour( labeledImage )
% Input:
%   labeledImage: labeled or BW image.
% Output:
%   contour: labeled / BW contour of the input image.
% Example:
%     I = imread('coins.png'); I = I(20:75, 30:85);
%     BW = im2bw(I);
%     contour = CalcContour(BW);
%     f = figure; 
%     subplot(3,1,1); imagesc(I); title('Intensity Image'); colormap gray;
%     subplot(3,1,2); imagesc(BW); title('Input binary image');
%     subplot(3,1,3); imagesc(contour); title('Output contour image');
%     set(f, 'Position', get(0, 'Screensize'));
% See also: im2bw; bwlabel
%
% T. Gilad, 2017
%%
    contour = abs(labeledImage - imerode(labeledImage,ones(3)));
end

