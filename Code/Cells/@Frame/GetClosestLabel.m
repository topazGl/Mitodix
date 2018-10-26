function [ label ] = GetClosestLabel( labeled, com )
% Get the label of the closest pixel to com.
%   label = GetClosestLabel( labeled, com )
% Inputs:
%   labeled = Labeled image.
%   com = 2x1 (X, Y) point.
% Outputs:
%   label = output label. 0 if there are no labels in the input image.
    [~, ind] = bwdist(labeled~=0);
    if(max(ind(:)))
        ind = ind(com(2), com(1));
        [cY, cX] = ind2sub(size(labeled), ind);
        label = labeled(cY, cX);
    else
        label = 0;
    end
end

