function cell = GetCell(self, cellLabel)
% Get Cell object of cellLabel
%     GetCell(self, cellLabel)
% Inputs:
%   self = Frame object.
%   cellLabel = cell label.
    if(self.cells.isKey(cellLabel))
        cell = self.cells(cellLabel);
    else
        [isOnBorder, self.PaddedImage, self.PaddedLabeled] = self.EstimateBorderCell(cellLabel);
        if(isOnBorder)
            self.setCellData(cellLabel);
        end
        cell = CellInstance(self, cellLabel);
        self.cells(cellLabel) = cell;
    end
end

