function [motherLabeles, foundFlags, R] = GetPossibleMothers(self, daughters, R, anchorMothers)
% Get possible mothers in frame, that are R radius away from the middle
% point of the input candidate daughters.
%   motherLabeles = GetPossibleMothers(self, daughtersMidPoints, R, anchorMothers)
% Inputs:
%   self = Frame object (mother frame).
%   daughtersMidPoints = Nx2 matrix of neigbours mid-points.
%   R = (optional) search radius [pixels].
% Output:
%   motherLabeles = Nx1 cell array. Each cell is a vector of possible
%                   mother labeles.
%   foundFlags = Nx1 array of flags. n each location: 1 - if mothers were
%                found, else - 0.
%   R = search radius.
    if(isempty(daughters))
        motherLabeles = [];
        foundFlags = 0;
        return;
    end
    
    N = size(daughters, 1);
    motherLabeles = cell(N, 1);
    foundFlags = zeros(N, 1);
    
    % Calculate default search radius
    if(~exist('R', 'var') || isempty(R) || (R <= 0))
        R = self.GetDriftRadius();
    end
    
    closeLabeles = [];
    if(exist('anchorMothers', 'var') && ~isempty(anchorMothers))
        % mothers must be near the anchors:
        for m=1:length(anchorMothers)
            closeLabeles = [closeLabeles; self.GetCellLabelsInR(anchorMothers(m), 2*R, 0)];
        end
    end
        
    for n=1:N
        % Get cells whitin the search radius:
        mLabeles = self.GetCellLabelsInR(daughters(n, :), R, 0);
        
        if(~isempty(closeLabeles))
            % mothers must be near the anchors:
            tmpLabeles = [];
            for m=1:length(mLabeles)
                if(any(closeLabeles == mLabeles(m)))
                    tmpLabeles = [tmpLabeles, mLabeles(m)];
                end
            end
            mLabeles = tmpLabeles;
        end
        
        if(isempty(mLabeles))
            return;
        end
        
        %Only consider cells that are in the right stage of the
        %cell cycle (ignor this condition in the first frames):
        ages = self.GetAges(mLabeles);
        cycle = self.Configure.CellCycle;
        if(~isempty(ages) && (self.Number >= cycle))
            idx = find(ages > 0.7*(23/24)*cycle);
            if(isempty(idx))
                return;
            end
            mLabeles = mLabeles(idx);
        end
        
        motherLabeles{n} = mLabeles;   
        foundFlags(n) = ~isempty(mLabeles);
        foundFlags = foundFlags>0;
    end

end

