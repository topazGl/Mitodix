function [ labels, scores ] = GetAnchorCells( self, refFrame, R, Rm )
    labels = [];
    
    minSize = floor(0.25*self.CellSize.Min);
    
    I1 = double(refFrame.Image);
    I2 = double(self.Image);
    
    % enhance low SNR cells:
    fg = I2(self.Labeled > 0);
    avg = mean(fg); %max(mean(fg) - std(fg), 0);
    L = length(self.labels);
    for l=1:L
        label = self.labels(l);
        M = mean(I2(self.Labeled == label));
        factor = avg/M;
        I2(self.Labeled == label) = factor*I2(self.Labeled == label);
    end
  
    fg = I1(refFrame.Labeled > 0);
    avg = mean(fg);%max(mean(fg) - std(fg), 0);
    L = length(refFrame.labels);
    for l=1:L
        label = refFrame.labels(l);
        M = mean(I1(refFrame.Labeled == label));
        factor = avg/M;
        I1(refFrame.Labeled == label) = factor*I1(refFrame.Labeled == label);
    end
    
    [locations_d, scores_d] = Frame.AttentionAnchors(I1, I2, minSize);
    
    % Also get hints from the opposite direction:
    [locations_m, scores_m] = Frame.AttentionAnchors(I2, I1, minSize);
    
    locations = [locations_d; locations_m];
    scores = [scores_d; scores_m];
    if(isempty(locations))
        return;
    end
    
    idx = 1:length(scores);
    if(exist('R', 'var') && ~isempty(R))
        if(size(locations, 1) == 1)
            return;
        end
        
        [~, D] = knnsearch(locations, locations, 'K', 2);
        D = D(:,2);

        idx = find(D <= R);
    end
    
    idx_m = idx(idx > length(scores_d));
    scores_m = scores(idx_m);
    
    idx_d = idx(idx <= length(scores_d));
    labels = self.GetClosestCellLabel(locations(idx_d, :));
    scores = scores(idx_d);

    if(~exist('Rm', 'var') || isempty(Rm))
        Rm = refFrame.GetDriftRadius();
    end
    
    for j=1:length(idx_m)
        labels_m = self.GetCellLabelsInR(locations(idx_m(j),:), Rm);
        labels = [labels, labels_m];
        scores = [scores; scores_m(j)*ones(length(labels_m), 1)];
    end
    
    [labels, idx] = unique(labels);
    scores = scores(idx);   
end

