function [labels, dists] = GetCellLabelsInR(self, ref, R, k)
% Find the k-nearest cell labels in R around the ref.
%   [labels, dists] = GetCellLabelsInR(self, ref, R, k)
% Inputs:
%   self = Frame object.
%   ref = a cell label in this frame (scalar), or a coordinate (2x1).
%   R = search radius [pixels].
%   k = k-nearest. If k=0 get all.
% Ouputs:
%   labels = output cell labeles.
%   dists = distances of output cells from ref.
    labels = [];
    dists = [];
    try
        if (length(ref) == 2)
            % ref is coordinates
            midPoint = ref;
        else
            % ref is label
            midPoint = self.GetCellData(ref);
        end

        if((nargin > 3) && (~k))
            k = length(self.labels);
        end

        % Check that input R is valid:
        if ((nargin > 3) && isnan(R))
            error(['Search radius R must be a positive number (input was ', num2str(R),').']);
        end

        if ((nargin < 3)  || (R <= 0))
            if(nargin > 3)
                % k-nearest:
                idx = knnsearch(self.centers, midPoint, k);
            else
                idx = knnsearch(self.centers, midPoint);
            end
        else
            % Only consider cells in radius R:
            R = round(R);
            idx = rangesearch(self.centers, midPoint, R);
            idx = idx{:};
            if (nargin > 3 && (length(idx) > k) && (k>0))
                % k-nearest:
                idx = idx(1:k);
            end
        end

        if (isempty(idx))
            return;
        end

        labels = self.labels(idx);
        
        if(nargout > 1)
            COMs = self.centers(idx,:);
            dists =  sqrt((COMs(:,1)-midPoint(1)).^2 + (COMs(:,2)-midPoint(2)).^2);
            if(length(dists) > 1)
                % sort according to distance:
                [dists, I] = sort(dists);
                labels = labels(I);
            end
        end
    catch EX
        self.Log.mlog = {LogType.ERROR ...
                        ,mfilename('class') ...
                        ,[self.Log.Me,' Problem finding close cells in frame ', num2str(self.Number), ': ', EX.message]};
        throw(EX);
    end
end

