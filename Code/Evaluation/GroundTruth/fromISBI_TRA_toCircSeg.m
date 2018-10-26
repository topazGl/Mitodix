function fromISBI_TRA_toCircSeg( dataName, outputPath, radi )
    if(~exist('radi', 'var') || isempty(radi))
        radi = 50;
    end
    
    if(~exist(outputPath, 'dir'))
        mkdir(outputPath);
    end
    mtdx = Mitodix(0, dataName);
    [labeldList, startFrame, endFrame] = mtdx.GetFrameImages(0, 0, 1);
    N = endFrame-startFrame+1;
    lblPath = mtdx.InputSegPath;
    for n = 1:N
    	t = n + startFrame - 1;
        labeled = uint16(imread(fullfile(lblPath, labeldList{n})));
        s = regionprops(labeled,'Centroid');
        centers = round(cat(1, s.Centroid));
        centers = centers(centers(:,1)>0, :);
        [Y, X] = size(labeled);
        [x, y] = meshgrid(1:X, 1:Y);
        new_seg = zeros(Y, X);
        for c=1:size(centers, 1)
            diskPixels = (y - centers(c,2)).^2 + (x - centers(c,1)).^2 <= radi.^2;
            new_seg(diskPixels) = labeled(centers(c,2), centers(c,1));
        end
        imwrite(new_seg, fullfile(outputPath, labeldList{n}));
    end
end

