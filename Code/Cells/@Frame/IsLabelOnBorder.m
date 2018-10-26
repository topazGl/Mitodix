function isOnBorder = IsLabelOnBorder(self,cellLabel)
% Check if input label is a cell on the border of the frame (8-connectivity).
%   isOnBorder = IsLabelOnBorder(self,cellLabel)
% Inputs:
%   self = Frame object.
%   cellLabel = cell label.
% Ouputs:
%   isOnBorder = cell on border flag.
    bw = (self.Labeled == cellLabel);
    filterBorderCell = imclearborder(bw, 4);
    isOnBorder = not(sum(filterBorderCell(:)));
end

