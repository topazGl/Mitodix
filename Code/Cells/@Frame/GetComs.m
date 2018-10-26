function coms = GetComs(self, margin)
% Get inner cell centers
%     coms = GetComs(self, margin)
% Inputs:
%   self = Frame object.
%   margin = number of margin pixels.
    coms = self.CellsData.Centers;
    if(margin == 0)
        return;
    end
    
    if(margin < 1)
        lbls = self.CellsData.Labels;
        innerLabels = self.Labeled;			
		innerLabels(imclearborder(innerLabels ~= 0) == 0) = 0;
        idx = [];
        for ii=1:length(lbls)
            if(any(innerLabels(:) == lbls(ii)))
                idx = [idx; ii];
            end
        end
    else
		Y = size(self.Image, 1);
		X = size(self.Image, 2);
        idx = ((coms(:,1) > margin) & (coms(:,2) > margin) & (coms(:,1) < X-margin) & (coms(:,2) < Y-margin));
    end
    coms = coms(idx,:);
end

