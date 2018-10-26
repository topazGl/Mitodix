function DisksSeg( dataName, outputPath, radi )
    if(~exist('radi', 'var') || isempty(radi))
        ellipse_flag = 1;
	else
		ellipse_flag = 0;
    end
    
    mtdx = Mitodix(0, dataName);
    mtdx.ReadFrames();
    endFrame = mtdx.LastFrameNum;
    startFrame = mtdx.FirstFrameNum;
    N = endFrame-startFrame+1;
    lblPath = mtdx.InputSegPath;
    
	if(~exist('outputPath', 'var') || isempty(outputPath))
		parts = strsplit(lblPath, filesep);
        outputPath = fullfile(parts{1:end-1}, 'DiskSeg');
	end
	
	if(~exist(outputPath, 'dir'))
        mkdir(outputPath);
    end
	
    [Y, X] = size(mtdx.LastFrame.Labeled);
    [x, y] = meshgrid(1:X, 1:Y);
    if(ellipse_flag)
        se = strel('disk',1);
        x0 = floor(X/2) + 1;
        y0 = floor(Y/2) + 1;
        RA = imref2d([Y, X]);
    end
    
    for n = 1:N
    	t = n + startFrame - 1;
        frame = mtdx.GetFrame(t);
        CellsData = frame.CellsData;
		
        new_seg = uint16(zeros(Y, X));
        for c=1:length(CellsData.Labels)
			label = CellsData.Labels(c);
            center_XY = CellsData.Centers(c, :);
			if(~ellipse_flag)
				diskPixels = (y - center_XY(2)).^2 + (x - center_XY(1)).^2 <= radi^2;
			else
				l = CellsData.Lengths(c)/2;
				w = CellsData.Widths(c)/2;
                psi = CellsData.Orientations(c);
				ellipse = ((y-y0)/l).^2 + ((x-x0)/w).^2 <= 1;
                ellipse = imrotate(ellipse, 90+psi, 'bicubic', 'crop');
                ellipse = imclose((ellipse ~= 0), se);
                ellipse = imtranslate(ellipse, RA, center_XY - [x0, y0], 'cubic');
                diskPixels = imclose((ellipse ~= 0), se);
			end
            new_seg(diskPixels) = label;
        end
        imwrite(new_seg, fullfile(outputPath, frame.LabeledName));
    end
	sprintf('Outputs saved to: %s', outputPath);
end

