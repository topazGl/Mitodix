function ages = GetAges(self, lbls)
% Get ages of input cells.
% Inputs:
%   self = Frame object.
%   lbls = array of input cell labeles.
% Output:
%   ages = array of cell ages [frames].
    if(~isempty(self.ages))
        ages = self.ages(ismember(self.labels, lbls));
    else
        ages = -1;
    end
end

