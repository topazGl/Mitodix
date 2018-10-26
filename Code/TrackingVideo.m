%%%
%	This function generates the tracking videos presented in:
%	Gilad et al, Fully Unsupervised Symmetry-Based Mitosis Detection in Time-Lapse Cell Microscopy, 2018.
%   (c) T.Gilad, 2018.
%%%
% Create an output video, mark tracks
function TrackingVideo(inDirIntensity, inDirTracking, expr, track_file_names, outDir, tracks, cellFrameRate, FOI_E, only_divisions_flag)
    %% Annotate frame by frame:
    vidWriter = VideoWriter.empty;
    if(~exist('w', 'var') || isempty(FOI_E))
        FOI_E = 0;
    end
    
    if(~exist('only_divisions_flag', 'var') || isempty(only_divisions_flag))
        only_divisions_flag = 0;
    end
    
    % Output folder:
    if(~exist(outDir, 'dir'))
        mkdir(outDir);
    end
    
    vidFolder = fullfile(outDir, 'Visualize');
    if(~exist(vidFolder, 'dir'))
        mkdir(vidFolder);
    end

    ext = {strfind(expr,'.')};
    ext = expr(ext{end}:end);
    fileName = dir(fullfile(inDirIntensity, ['*',ext]));
    [~, intensity_file_names] = Frame.GetFileNameToken({fileName.name}, expr);

    % choose a diferent color for each candidate:
    colorsNames = {'cyan', 'magenta', 'yellow', 'white'};
    colorsRGB = 255*[0, 1, 1; 1, 0, 1; 1, 1, 0; 1, 1, 1];
    colorsN = length(colorsNames);
    N = length(intensity_file_names);
    % Write frames into video:
    for n = 1:N
        I = double(imread(fullfile(inDirIntensity, intensity_file_names{n})));
        L = imread(fullfile(inDirTracking, track_file_names{n}));
        bgcolor = mean(I(L(:)==0));
        I = I - bgcolor;
        
        s = regionprops(L, 'Centroid');
        centers = cat(1, s.Centroid);
        labels = sort(unique(L(:)));
        labels = labels(labels>0);
        centers = centers(labels, :);
%         centers = centers(~isnan(centers(:, 1)), :);

        % read input frame image (gray scale):
        I = (2^8-1) * I/max(I(:));
        I = uint8(I);
        I = cat(3,I,I,I);
        fNum = n-1;
        for ll=1:length(labels)
            color_ind = 1;
            
            %check if new:
            cell_track = (tracks(:,1)==labels(ll));
            birth_track = (tracks(:,2)==fNum);
            if ((n>1) && any(cell_track & birth_track))
                if(any(cell_track & birth_track & tracks(:,4)~=0))
                    color_ind = 2;
                else
                    color_ind = 4;
                end
            else
                cell_as_parent = (tracks(:,4)==labels(ll));
                parent_track = (tracks(:,2)==fNum+1);
                if((n<N) && any(cell_as_parent & parent_track))
                    color_ind = 3;
                end
            end
        
            if((color_ind==1 || color_ind==4) && only_divisions_flag)
                continue;
            end
            
            color = colorsRGB(color_ind, :);

            %mark contour:
            idxs = CalcContour(L == labels(ll))~=0;
            IR = I(:,:,1); IG = I(:,:,2); IB = I(:,:,3);
            IR(idxs) = color(1); IG(idxs) = color(2); IB(idxs) = color(3);
            I = cat(3,IR,IG,IB);

            %mark COM:
            I = insertMarker(I, centers(ll,:) , 'x', ...
                                'Color', colorsNames{color_ind});
            
            txt = num2str(labels(ll));
            font_size=20;
            if( color_ind == 2 )
                track_ind = find(tracks(:,1)==labels(ll), 1, 'first');
                parent = tracks(track_ind, 4);
                if(parent)
                    if(only_divisions_flag)
                        txt = num2str(parent);
                    else
                        txt = [txt,  '_', num2str(parent)];
                        font_size=16;
                    end
                end
            end
            % Mark serial number:
            textLocation = int32([centers(ll,1), centers(ll,2)-5]);
            I = insertText(I,textLocation,txt,'FontSize',font_size,'TextColor',colorsNames{color_ind},'BoxOpacity',0);
        end
        % Add frame into video:
        
        Iclip = I(FOI_E+1:end-FOI_E,FOI_E+1:end-FOI_E, :);
        imwrite(Iclip, fullfile(vidFolder, ['Vis_', intensity_file_names{n}]));
    end

    try
        frameRate = round(20/cellFrameRate);

        % Create a video writer object:
        vidWriter = VideoWriter(fullfile(vidFolder, 'Tracking_Results.avi'));
        vidWriter.FrameRate = 10;
        repeat = round(vidWriter.FrameRate/frameRate);

        % Open the video writer object:
        open(vidWriter);

        expr = ['Vis_', expr];
        fileList = dir(fullfile(vidFolder, 'Vis_*'));                
        fileNames = {fileList.name};
        [t, fileNames] = Frame.GetFileNameToken( fileNames, expr );
        [t,sortID] = sort(t);
        fnames = fileNames(sortID);

        for kk=1:length(fnames)
            Iclip = imread(fullfile(vidFolder,fnames{kk}));
            for ii=1:repeat
                writeVideo(vidWriter,insertText(Iclip, [0 0; 0 0], ['t',num2str(t(kk))], 'FontSize', 24));
            end
        end
        close(vidWriter);
    catch EX
        if(~isempty(vidWriter))
            close(vidWriter);
        end
        rethrow(EX);
    end
end
