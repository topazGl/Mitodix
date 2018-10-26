function setCellData(self, cellLabel)
% Recalculate cell data.
%   setCellData(self, cellLabel)
% Inputs:
%   self = Frame object.
%   cellLabel = cell label.
    idx = find(self.labels == cellLabel, 1);
    if(idx <= 0)
        return;
    end
    msk = (self.PaddedLabeled == cellLabel);
    s = regionprops(msk, 'Centroid','MajorAxisLength', ...
        'MinorAxisLength', 'Area', 'Orientation');
    self.lengths(idx) = cat(1, s(1).MajorAxisLength);
    self.widths(idx) = cat(1, s(1).MinorAxisLength);
    self.sizes(idx) = cat(1, s(1).Area);
    self.centers(idx,:) = floor(cat(1, s(1).Centroid)) - [self.W, self.W];
    self.orientations(idx) = cat(1, s(1).Orientation);
end

