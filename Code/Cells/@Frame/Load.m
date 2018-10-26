function Load(self, sobj)
% Load struct data
%   Load(self, sobj)
% Inputs:
%   self = Frame object.
%   sobj = Frame struct.
    self.Number = sobj.Number;                  
    self.Name = sobj.Name;
    self.LabeledName = sobj.LabeledName;
    self.Image = sobj.Image;                      
    self.Labeled = sobj.Labeled;                  
    self.labels = sobj.CellsData.Labels;                 
    self.centers = sobj.CellsData.Centers;                  
    self.ages = sobj.CellsData.Ages;                      
    self.neighbors = sobj.CellsData.Neighbors;                
    self.neighborsDists = sobj.CellsData.NeighborsDists;   
    self.neighborsMidPoints = sobj.CellsData.NeighborsMidPoints; 
    self.lengths = sobj.CellsData.Lengths;                   
    self.widths = sobj.CellsData.Widths;                    
    self.sizes = sobj.CellsData.Sizes;
    self.orientations = sobj.CellsData.Orientations;
    self.Size = sobj.Size;                     
    self.CellLength = sobj.CellLength;                 
    self.CellWidth = sobj.CellWidth;                  
    self.CellSize = sobj.CellSize;                  
    self.CellDrift = sobj.CellDrift;                     
    self.DistBetweenNeigbors = sobj.DistBetweenNeigbors;
%     self.PaddedImage = sobj.PaddedImage;
%     self.PaddedLabeled = sobj.PaddedLabeled;
%     self.W = sobj.W;
%    self.cells = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    self.SaturationFlag = sobj.SaturationFlag;
	self.Z_dim = sobj.Z_dim;
end

