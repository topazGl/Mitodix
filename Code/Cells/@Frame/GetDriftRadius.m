function R = GetDriftRadius( self )
% Get search radius based on cell drift statistics.
%   R = GetDriftRadius( self )
% Input:
%   self = Frame object.
% Output:
%   R = search radius [pixels].
    if(self.NumberOfCells < 50)
        % No statistics
        R = inf;
        return;
    end
    
    if(~isempty(self.CellDrift))
        mu = round(self.CellDrift.Mean);
        s = round(self.CellDrift.STD);
    else
        mu = 0;
        s = 0;
    end
    l = 0.5*max(self.CellLength.Min, self.CellLength.STD);
    
    if(s <= mu)
        R = ceil(mu + 2*s + l);% 3*s + 2*l);
    else
        if(l > (mu+s))
            R = ceil(mu + s + 2*l);
        else % cell length variation neglected:
            R = ceil(mu + s);
        end
    end
end

