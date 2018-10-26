function [center, length, width, size, orientation] = GetCellData(self, cellLabel)
% Get the data for the input cell label.
%   [center, length, width, size, orientation] = GetCellData(self, cellLabel)
% Inputs:
%   self = Frame object.
%   cellLabel = cell label.
% Outputs:
%   center = cell COM
%   length = cell length [pixels].
%   width = cell width [pixels].
%   size = cell size [pixels].
%   orientation = cell major axis orientation [deg].
    center = [];
    length = 0;
    width = 0;
    size = 0;
    orientation = 0;
    if(~isempty(self.labels) && ~isempty(cellLabel))
        idx = find(self.labels == cellLabel,1);
        if(~isempty(idx) && (idx > 0))
            center = double(self.centers(idx, :));
            length = double(self.lengths(idx));
            width = double(self.widths(idx));
            size = double(self.sizes(idx));
            orientation = double(self.orientations(idx));
        end
    end
end

