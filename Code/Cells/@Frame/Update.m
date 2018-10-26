function Update(self, varargin)
% Update frame data.
%   Update(self, varargin)
% Inputs:
%   self = Frame object.

    % Update cell ages only:
    if(nargin == 1)
        self.SetAges(varargin{1});
        return;
    end

    % Update cell labels and ages:
    self.cells = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    inLabels = varargin{1};
    newCOMs = varargin{2};
    [newLabels, idx] = GetClosestCellLabel(self, newCOMs, inLabels);
    newLabeled = zeros(self.Size);
    newCenters = zeros(1, length(idx));
    for n = 1:length(idx)
        newLabeled(self.labeled == self.labels(idx(n))) = newLabels(idx(n));
        newCenters(n) = newCOMs(idx(n));
    end
    newAges = varargin{3}(inLabels == newLabels);

    newNeigbours = zeros(length(self.neigbours), 2);
    for n=1:length(self.neigbours)
        newNeigbours(n,:) = [newLabeled(self.labeled == newNeigbours(n,1)), newLabeled(self.labeled == newNeigbours(n,1))];
    end
    self.labeles = newLabels;
    self.labeled = newLabeled;
    self.centers = newCenters;
    self.neigbours = newNeigbours;
    self.SetAges(newAges);
end

