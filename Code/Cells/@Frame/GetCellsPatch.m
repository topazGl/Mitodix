function [ patch, weights, mask, cellsImages, coms, d_deg ] = GetCellsPatch( self, labels, pSize, alignFlag, divideFlag )
% Get patch around midPoint containing only cells of requested 2 labels where the cells are one above the other.
%     [ patch, mask, weights ] = GetCellsPatch( self, labels )
% Inputs:
%   self = Frame object.
%   labels = 2x1 cells labeles in frame.
%   pSize = optional. Patch size [pSize, pSize], pSize = 2*w+1.
%   alignFlag = Adjust cells alignment according to symmetry axis (Default = 0).
%   divideFlag = Create additional patch images for each cell (Default = 0).
%
% Ouputs:
%   patch = [2*w+1, 2*w+1] patch image of requested labels.
%   weights = [2*w+1, 2*w+1] values 0-1 of probability of the pixel to
%             belong to one of the cells based on segmentation.
%   mask = [2*w+1, 2*w+1] binary mask of requested cells ('1' = cell pixel).
%   cellsImages = 2x1 cell array of cells patch images (if divideFlag = 0:
%                 value is empty).
%   coms = 2x2 cells COMs in the frame.
%   d_deg = rotation angle
    cellsImages = [];
    
    if(isempty(labels) || length(labels) ~= 2)
        error('GetCellsPatch Must get 2 labels!');
    end
    
    if(~exist('alignFlag', 'var') || isempty(alignFlag))
        alignFlag = 0;
    end
    
    if(~exist('divideFlag', 'var') || isempty(divideFlag))
        divideFlag = 0;
    end
    
    d_deg = 0;
    coms = ones(2,3);
    [com, l1] = self.GetCellData(labels(1));
    coms(1,1:length(com)) = com;
    [com, l2] = self.GetCellData(labels(2));
    coms(2,1:length(com)) = com;
    l = max(l1, l2);
    midPoint = 0.5*(coms(1,:) + coms(2,:));
    
    if(exist('pSize', 'var') && ~isempty(pSize))
        if(mod(pSize, 2) == 0)
            w = pSize/2;
        else
            w = (pSize - 1)/2;
        end
    else
        d = norm(com1 - com2);
        w = round(0.5*d + l);
    end
    
	[lbld, I] = self.GetPadded(w, 0, round(midPoint(3)));
    midPoint = round(midPoint(1:2)) + [w, w];
    
    subLabeled = double(lbld(midPoint(2)-w:midPoint(2)+w, midPoint(1)-w:midPoint(1)+w));
    patch = double(double(I(midPoint(2)-w:midPoint(2)+w, midPoint(1)-w:midPoint(1)+w)));

    mask1 = (subLabeled == labels(1));
    mask2 = (subLabeled == labels(2));
    
    mask = mask1 | mask2;
    subLabeled(~mask) = 0;
%     s = regionprops(subLabeled,'Centroid');
%     subCenters = cat(1, s.Centroid);
%     subCenters = subCenters(~isnan(subCenters(:,1)) & ~isnan(subCenters(:,2)), :);
    
    s = regionprops(mask1,'Centroid', 'Area');
    subCenters1 = cat(1, s.Centroid);
    subAreas = cat(1, s.Area);
    [~, max_ind] = max(subAreas);
    subCenters = subCenters1(max_ind, :);
    
    s = regionprops(mask2,'Centroid', 'Area');
    subCenters2 = cat(1, s.Centroid);
    subAreas = cat(1, s.Area);
    [~, max_ind] = max(subAreas);
    subCenters = [subCenters; subCenters2(max_ind, :)];
    
     if(size(subCenters, 1) < 2)
        warning(['GetCellsPatch: In frame #', num2str(self.Number), ' patch size of cells ', num2str(labels(1)), ' and ', num2str(labels(2)), ' is too small!']);
        patch = [];
        mask = [];
        weights = [];
        d_deg = 0;
        return;
     end
    
    % initial registration (segmentation based):
    middle = mean(subCenters, 1);
    dt = [w+1, w+1] - middle;
    d_deg = Frame.CalcOrientation(subCenters(1,:), subCenters(2,:));
    [patch, subLabeled, mask, weights] = DaughtersImTransform(dt, d_deg, patch, subLabeled, labels);
    
    % align according to symmetry axis:
    if(alignFlag)
        x0 = [0, 0, 0];
        x_limit = [0.25*l, 0.25*l, 15];
        options_fminsearch = optimset('TolX',0.1, 'TolFun', 0.0001, 'Display','off');
        x = fmincon(@(x) DaughtersDissimilarity(patch, subLabeled, labels, x), x0, [], [], [], [], -x_limit, x_limit,[], options_fminsearch);
        
        dt = x(1:2)/2;
        deg_corr = x(3)/2;
        [patch, subLabeled, mask, weights] = DaughtersImTransform(dt, deg_corr, patch, subLabeled, labels);
    else
        deg_corr = 0;
    end
    d_deg = d_deg + deg_corr;
    
    % create cells patches:
    if(divideFlag)
        cellsImages = cell(2, 1);
    
        pLabeled = padarray(subLabeled,[w,w], 'both');
        pIm = padarray(patch,[w,w], 'both');
        
        s = regionprops(pLabeled,'Centroid');
        subCenters = cat(1, s.Centroid);
        subCenters = subCenters(~isnan(subCenters(:,1)) & ~isnan(subCenters(:,2)), :);
        subCenters = round(subCenters);
        
        if(size(subCenters, 1) ~= 2)
            if(size(subCenters, 1) < 2)
            	warning(['GetCellsPatch: In frame #', num2str(self.Number), ' problem dividing cells ', num2str(labels(1)), ' and ', num2str(labels(2)), ', less than 2 cells within sub-patch!']);
            else
                warning(['GetCellsPatch: In frame #', num2str(self.Number), ' problem dividing cells ', num2str(labels(1)), ' and ', num2str(labels(2)), ', more than 2 cells within sub-patch!']);
            end
            cellsImages = [];
            d_deg = 0;
            return;
        end
        
        for jj=1:2
            cellsImages{jj}.Im = double(pIm(subCenters(jj, 2)-w:subCenters(jj, 2)+w, subCenters(jj, 1)-w:subCenters(jj, 1)+w));
            cellsImages{jj}.MaskIm = double(pLabeled(subCenters(jj, 2)-w:subCenters(jj, 2)+w, subCenters(jj, 1)-w:subCenters(jj, 1)+w));
            cellsImages{jj}.MaskIm = (cellsImages{jj}.MaskIm == labels(jj));

            [~,~,~,pW] = Frame.TransformLabeledImage(cellsImages{jj}.MaskIm, 1);
            maxMask = max(pW(:));
            if(maxMask ~= 0)
                pW = double(pW)/max(pW(:));
            end
            cellsImages{jj}.Weights = pW;
            cellsImages{jj}.Masked = cellsImages{jj}.Im .* pW;
            
            if(jj == 2)
                cellsImages{jj}.Im = flipud(cellsImages{jj}.Im);
                cellsImages{jj}.MaskIm = flipud(cellsImages{jj}.MaskIm);
                cellsImages{jj}.Weights = flipud(cellsImages{jj}.Weights);
                cellsImages{jj}.Masked = flipud(cellsImages{jj}.Masked);
            end
        end
    end
    coms = coms(:, 1:2);
end

function [patch, labeled, mask, weights] = DaughtersImTransform(dt, d_deg, patch, labeled, labels)
    if(~exist('labels', 'var') || isempty(labels))
        labels = unique(labeled(:));
        labels = labels(labels ~= 0);
    end
    
    patch = RigidTransform(patch, dt, d_deg);
    mask1 = RigidTransform(labeled == labels(1), dt, d_deg, 1);
    mask2 = RigidTransform(labeled == labels(2), dt, d_deg, 1);
    labeled = zeros(size(labeled));
    labeled(mask1) = labels(1);
    labeled(mask2) = labels(2);
    labeled(mask1 & mask2) = 0;
    
    mask = labeled ~= 0;
    
    [~,~,~,weights] = Frame.TransformLabeledImage(mask, 1);
    maxMask = max(weights(:));
    if(maxMask ~= 0)
        weights = double(weights)/max(weights(:));
    end
end

function dissimilarity = DaughtersDissimilarity(patch, subLabeled, labels, x)
    dt = x(1:2);
    d_deg = x(3);
    fixedIm = flipud(patch);
    fixedMask = flipud(subLabeled ~= 0);
    
    [movingIm, ~, movingMask] = DaughtersImTransform(dt, d_deg, patch, subLabeled, labels);
    
    [~,~,~,weights] = Frame.TransformLabeledImage(movingMask | fixedMask, 1);

    errVec = fixedIm(:) - movingIm(:);
    dissimilarity = errVec'*(weights(:).*errVec);
%     dissimilarity = -wcorr2( fixedIm, movingIm, weights );
end
