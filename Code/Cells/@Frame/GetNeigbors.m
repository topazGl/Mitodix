function [neighborsLabels, neighborsMidPoints, neighborsDists] = GetNeigbors(self, R, cellLabels)
% Get neigbor cell data of input cells that are whitin an R radius from each other.
% Inputs:
%   self = Frame object.
%   R = search radius [pixels].
%   cellLabels = array of cell labels.
% Outputs:
% 	neighborsLabels = Mx2 matrix of all neigboring cell pairs (cell labels).
%	neighborsDists =  Mx1 vector of all neigboring cell pairs distances.
%	neighborsMidPoints =  Mx2 matrix of all neigboring cell mid points.
    neighborsLabels = [];
    neighborsMidPoints = [];
    neighborsDists = [];
    if(nargin < 3)
        inLabels = 0;
    else
        inLabels = 1;
    end

    try
        % Check that input R is valid:
        if ((R <= 0) || isnan(R))
            error(['Search radius R must be a positive number (input was ', num2str(R),').']);
        end

        idx = (self.neighborsDists <= R)';
        closeNeigbors = self.neighbors(idx,:);
        closeNeigborsD = self.neighbors(idx);
        closeNeigborsM = self.neighborsMidPoints(idx, :);
        
        N = size(closeNeigbors,1);
        for n = 1:N
            pair = closeNeigbors(n,:);

            % consider only neigbors of input cells:
            if (inLabels && not(any((cellLabels == pair(1))|(cellLabels == pair(2)))))
                continue;
            end 
            neighborsLabels = [neighborsLabels; pair];
            neighborsMidPoints = [neighborsMidPoints; closeNeigborsM(n,:)];
            neighborsDists = [neighborsDists, closeNeigborsD(n)];
        end

   catch EX
        self.Log.mlog = {LogType.ERROR ...
                        ,mfilename('class') ...
                        ,[self.Log.Me,' Problem finding negibors in frame ', num2str(self.Number), ': ', EX.message]};
        throw(EX);
    end
end

