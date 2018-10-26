function SeparateCellsSegBatch(mtdx, startSegFrame, endSegFrame)
% Perform segmentation re-calculation on an images batch
%   SeparateCellsSegBatch(mtdx, startSegFrame, endSegFrame)
% Inputs:
%   mtdx = Mitodix object.
%   startSegFrame = start batch frame number.
%   endSegFrame = end batch frame number.

    % Get input segmentation:
    [labeldList, startFrame, endFrame] = mtdx.GetFrameImages(startSegFrame, endSegFrame, 1);
    N = endFrame-startFrame+1;
    if(~N)
        return;
    end

    % Output directory:
    newSegFolder = fullfile(mtdx.OutputPath,'newSeg');
    if(~exist(newSegFolder, 'dir'))
        mkdir(newSegFolder);
    end

    % Save new segmentation:
    for n = 1:N
        seg = round(imread(fullfile(mtdx.InputSegPath, labeldList{n})));
        seg = seg - min(seg(:)); % back ground label must be zero.
        bw = Frame.SeparateCellsSeg(seg);
        imwrite(uint16(bw), fullfile(newSegFolder, labeldList{n}));
    end
end

