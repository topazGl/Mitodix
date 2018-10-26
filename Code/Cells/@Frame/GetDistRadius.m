function R = GetDistRadius( self )
% Get search radius based on cell neigbhors statistics.
%   R = GetDistRadius( self )
% Input:
%   self = Frame object.
% Output:
%   R = search radius [pixels].
    dists = self.neighborsDists(self.neighborsDists < self.DistBetweenNeigbors.Mean);
    dists = sort(dists);
    
    mu = mean(dists);
    s = std(dists);
    m = min(dists);
    med = dists(max(floor(length(dists)/2), 1));
    mu = min(mu, med);
    
    l = self.CellLength.Mean + 2*self.CellLength.STD;
    if(s/mu > 0.5)
        R = ceil(min(max([mu - s, m + 0.5*l, 2*m]), mu));
    else
        R = ceil(min(max([mu - 2*s, m + 0.5*l, 2*m]), mu));
    end
%     R = ceil(min(max(mu - 2*s, m + 0.5*l), mu));
    % R = ceil(min(max(mu - s, m + l), mu));
end

