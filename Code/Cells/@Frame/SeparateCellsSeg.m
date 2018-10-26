function bw = SeparateCellsSeg(seg)
% Re-calculate segmentation to separate connected cells
%   bw = SeparateCellsSeg(seg)
% Input:
%   seg = input segmentation image.
% Output:
%   bw = output binary image.

    bw1 = (seg ~= 0);
    labeled = (labelmatrix(bwconncomp(bw1)));

    % find connected cells:
    s = regionprops(labeled, 'Solidity');
    solidity = cat(1, s.Solidity);
    doubleCellsLabels = find(solidity < 0.85);
    doubleCellsIdx = [];
    for i=1:length(doubleCellsLabels)
        idx = find(labeled == doubleCellsLabels(i));
        doubleCellsIdx = [doubleCellsIdx; idx];
    end

    % separate using erode:
    doubleCells = zeros(size(bw1));
    doubleCells(doubleCellsIdx) = labeled(doubleCellsIdx);
    se = strel('disk',2);
    doubleCellsBW = imerode(doubleCells~=0, se);

    bw = bw1;
    bw(doubleCellsIdx) = doubleCellsBW(doubleCellsIdx);

    % check if any still connected:
    labeled = labelmatrix(bwconncomp(doubleCellsBW));
    s = regionprops(labeled, 'Solidity');
    solidity = cat(1, s.Solidity);
    doubleCellsLabels = find(solidity < 0.9);

    % separate connected cells:
    for i=1:length(doubleCellsLabels)
        b = zeros(size(bw));
        idx = (labeled == doubleCellsLabels(i));
        b(idx) = 1;
        D = -bwdist(~b);
        D(~b) = -Inf;
        L = watershed(D, 4);
        b(L == 0) = 0;
        b = imclose(b, se);
        bw(idx) = b(idx);
    end
end

