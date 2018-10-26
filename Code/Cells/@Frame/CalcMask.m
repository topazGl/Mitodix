function mask = CalcMask(BW)
% Calculate mask (round edges).
%   mask = CalcMask(BW)
% Input:
%   BW = binary input mask image.
% Ouput:
%   mask = binary output image (rounded edges).
    BW = imfill(BW,'holes');
    mask = imerode(BW, strel('disk', 1));
    mask = imdilate(mask, strel('disk', 2));
    if ~any(mask(:))
        mask = BW;
    end
end

