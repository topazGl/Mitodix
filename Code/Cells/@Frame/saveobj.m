function sobj = saveobj(self)
% Perform subclass save
%   sobj = saveobj(self)
% Input:
%   self = Frame object.
% Output:
%   sobj = Frame struct.
    sobj.Number = self.Number;                  
    sobj.Name = self.Name;
    sobj.LabeledName = self.LabeledName;
    sobj.Image = self.Image;                      
    sobj.Labeled = self.Labeled;                  
    sobj.CellsData = self.CellsData;                 
    sobj.Size = self.Size;                     
    sobj.CellLength = self.CellLength;                 
    sobj.CellWidth = self.CellWidth;                  
    sobj.CellSize = self.CellSize;                  
    sobj.CellDrift = self.CellDrift;                     
    sobj.DistBetweenNeigbors = self.DistBetweenNeigbors;
%     sobj.PaddedImage = self.PaddedImage;
%     sobj.PaddedLabeled = self.PaddedLabeled;
%     sobj.W = self.W;
    sobj.SaturationFlag = self.SaturationFlag;
	sobj.Z_dim = self.Z_dim;
end

