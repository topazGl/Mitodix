function [newCenters, newLabels] = GetNewerCellsCenters(self, maxPreviousLabelValue)
% Get cells with labeles bigger than the maximal in prev frame.
%   [newCenters, newLabels] = GetNewerCellsCenters(self, maxPreviousLabelValue)
% Inputs:
%   self = Frame object.
%   maxPreviousLabelValue = maximal label in the previous frame.
% Outputs:
%   newCenters = centers of all the cells in this frame that
%                have a label number bigger than maxPreviousLabelValue.
%   newLabels = The label values of those cells.
    idx = (self.labels > maxPreviousLabelValue);
    if(any(idx))
        newLabels = self.labels(idx);
        newCenters = self.centers(idx, :);
    else
        newLabels = [];
        newCenters = [];
    end
end

