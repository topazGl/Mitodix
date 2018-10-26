function SetCellDrift(self, prevFrameCenters)
% set cell drift according to minimal centers distances from 
% previous frame's centers.
%   updateFrame = SetCellDrift(self, prevFrameCenters)
% Inputs:
%   self = Frame object.
%   prevFrameCenters = 2D matrix of cells cenetrs in the previous frame.
    try
        [~,drifts] = knnsearch(self.centers,prevFrameCenters);
        self.CellDrift = Statistics(drifts);
    catch EX
        self.Log.mlog = {LogType.ERROR ...
                      ,mfilename('class') ...
                      ,[self.Log.Me,' Problem celculating CellDrift for frame ', num2str(self.frameNum),'.']};
        rethrow(EX);
    end
end

