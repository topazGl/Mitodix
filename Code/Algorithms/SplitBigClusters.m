function [ newLabeled, split ] = SplitBigClusters( labeled, std_factor )
    if(~exist('std_factor', 'var') || isempty(std_factor))
        std_factor = 3;
    end
    
    [ newLabeled, split ] = WaterShed( labeled );
    CC = bwconncomp(labeled~=0);
    big_labeled = labelmatrix(CC);
    s = regionprops(CC,'Area');
    areas = cat(1, s.Area);
    mean_area = mean(areas);
    sig_areas = std(areas);
    big_area = mean_area + std_factor*sig_areas;
    idx = find(areas <= big_area);
    if(~isempty(idx))
        for ii=1:length(idx)
            big_labeled(big_labeled == idx(ii)) = 0;
        end
    end
    [ ~, split_big ] = WaterShed( big_labeled );
    newLabeled(split_big) = 0;
    split = (split | split_big);
end

