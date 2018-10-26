function [labels, idx] = GetClosestCellLabel(self, newCOMs, newLabels)
% Get the label of the cell that is closest to newCOMs.
%   [labels, idx] = GetClosestCellLabel(self, newCOMs, newLabels)
% Inputs:
%   self = Frame object.
%   newCOMs = Nx2 array of cells COMs.
%   newLabels = (optional) Nx1 array of cells labeles.
% Outputs:
%   labeles = output labeles.
%   idx = labeles locations in array.
    if(isempty(self.centers) || isempty(newCOMs))
        labels = [];
        idx = [];
    else
        if (nargin > 2)
            idx = knnsearch(newCOMs, self.centers);
            labels = newLabels(idx);
        else
            idx = knnsearch(self.centers,newCOMs);
            labels = self.labels(idx);
        end
    end
end

