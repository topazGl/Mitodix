function [ GT_formatted ] = convertISBIGTformat( GT, manTrackPath, downscaleFactor, dt )
    if(ischar(GT))
        % read GT from path:
        GT = load(fullfile(GT, 'mitosisGT.mat'));
        GT = GT.mitosisGT;
    end
    if(~exist('downscaleFactor', 'var') || isempty(downscaleFactor))
        downscaleFactor = 1;
    end
    if(~exist('dt', 'var') || isempty(dt))
        dt = 0;
    end
    
    setCOM = 0;
    if(exist('manTrackPath', 'var') && ~isempty(manTrackPath))
        setCOM = 1;
        % Find all image files in input folder:
        fileList = dir(fullfile(manTrackPath,'*.tif'));
        [t, fileNames] = Frame.GetFileNameToken( {fileList.name}, 'man_track([0-9]+)' );
        if (isempty(t))
            [t, fileNames] = Frame.GetFileNameToken( {fileList.name}, 'mask([0-9]+)' );
        end
        if (isempty(t))
            [t, fileNames] = Frame.GetFileNameToken( {fileList.name}, 't([0-9]+)' );
        end
        
        [~,sortID] = sort(t);
        frameNamesList = fileNames(sortID);
        N = length(frameNamesList);

        % Check if input folder is not empty:
        if(~N)
            errMsg = ['There are no relevant files in frames input folder ', manTrackPath,'.'];
            error(errMsg);
        end
        
        parfor n=1:N
            [ ~, labeled ] = Frame.Read2D3D( fullfile(manTrackPath, fileNames{n}), fullfile(manTrackPath, fileNames{n}), [], [], downscaleFactor );
            s = regionprops(labeled, 'Centroid','Area');
            
            areas = cat(1, s.Area);
            labels{n} = find(areas > 0);
            centers{n} = round(cat(1, s(labels{n}).Centroid));
        end
    end
    
    GT_formatted = [];
    for d=1:size(GT,1)
        d1 = GT{d,1};
        d2 = GT{d,2};
        m = GT{d,3};
        cand_data.BirthFrame = d1.frameIndex;
        if(isfield(d1, 'center'))
            cand_data.MotherCOM = m.center;
            cand_data.Daughter1COM = d1.center;
            cand_data.Daughter2COM = d2.center;
        elseif(setCOM)
            idx = find(labels{d1.frameIndex-1-dt} == m.id, 1);
            cand_data.MotherCOM = centers{d1.frameIndex-1-dt}(idx(1), :);
            idx = find(labels{d1.frameIndex-dt} == d1.id, 1);
            cand_data.Daughter1COM = centers{d1.frameIndex-dt}(idx(1), :);
            idx = find(labels{d1.frameIndex-dt} == d2.id, 1);
            cand_data.Daughter2COM = centers{d1.frameIndex-dt}(idx(1), :);
        end
        if(isfield(d1, 'id'))
            cand_data.MotherLabel = m.id;
            cand_data.Daughter1Label = d1.id;
            cand_data.Daughter2Label = d2.id;
        end
        
        GT_formatted = [GT_formatted; cand_data];
    end
    GT_formatted = struct2table(GT_formatted);
end

