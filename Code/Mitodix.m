%%%
%	This is the main class for the Mitosis detection presented in:
%	Gilad et al, Fully Unsupervised Symmetry-Based Mitosis Detection in Time-Lapse Cell Microscopy, 2018.
% (c) T.Gilad, 2018.
%%%
classdef (Sealed) Mitodix < handle
% The main Mitosis detection class class.

%% Properties
    properties (Access = private)
        dataName = '';                  % Data set name (as in config file).
        frames;                         % Hashtable of Frame objects with
                                        % frame number as key.
        creationTime = '';              % creation time string of Mitodix
                                        % object.
                                        % Format: yyyy-mm-dd_HH-MM-SS-FFF.
        run_log = '';                   % specific run log file path.
        downscaleFactor = 1;            % resolution downscale factor.
        dt = 1;                         % frames gaps.
        frameNum = nan;                 % Used when DS is composed of several sites.
        cells_dists = [];               % vector of distances between nearest cells.
        cells_drifts = [];              % vector of cell drifts.
    end
    
    properties (SetAccess = private)
        Log;                            % Application logger.
        Configure;                      % Application configuration parameters.
    end
    
    properties
        DEBUG = 0;                      % Application Debug mode flag.
    end
    
    properties (Dependent)
        DataName;                       % data set name (as in config file).
        InputImPath;                    % input gray-level directory path.
        InputSegExpression;             % input segmentation file name expression.
        InputSegExtension;              % input segmentation files extension.
        InputSegPath;                   % input segmentation directory path.
        InputImExtension;               % input gray-level files extension.
        InputImExpression;              % input gray-level file name expression.
        OutputPath;                     % output directory name.  
        GTPath;                         % Ground truth *.mat file path.
        GroundTruth;                    % Ground Truth data (loaded from GTPath).
        MinCellLength;                  % min cell major axis (in pixels).
        FramesNumbers;                  % Array of all the frame numbers in the Mitodix collection.
        FirstFrameNum;                  % First frame number in Mitodix collection.
        LastFrameNum;                   % Last Frame object (most current).
        LastFrame;                      % Last Frame number (most current).
        Frames;                         % Array of all the Frame object in the Mitodix collection.
        Dists;                          % vector of distances between nearest cells.
        Drifts;                         % vector of cells drifts.
    end


%% Public methods
    methods
        % Mitodix Constructor
        function self = Mitodix(debugMode, dataName, downscaleFactor, dt, frameNum)
        % Mitodix Constructor
        % signature:
        %   mtdx = Mitodix(debugMode)
        % Input:
        %   debugMode 0/1 flag for debug outputs (Default = 0). Note: Debug 
        %   mode takes longer to execute.
        %   dataName: data set run name in config (optional). Default is
        %   according to chosen runName in config).
        %   downscaleFactor: 0<downscaleFactor<=1 resolution downscale
        %   (Default=1).
        %   dt: 1<=dt, frames gap (Default=1).
        %%
            if(~exist('downscaleFactor', 'var') || (isempty(downscaleFactor)))
                downscaleFactor = 1;
            end
            self.downscaleFactor = downscaleFactor;
            
            if(~exist('dt', 'var') || (isempty(dt)))
                dt = 1;
            end
            self.dt = dt;
            
            if(~exist('frameNum', 'var') || isempty(frameNum))
                frameNum = nan;
            end
            self.frameNum = floor(frameNum/self.dt);
            
            if(exist('dataName', 'var') && (~isempty(dataName)))
                % set chosen configuration:
                xmlFile = xmlread('Mitodix.config');
                allListItems=xmlFile.getElementsByTagName('root').item(0);
                thisListItem=allListItems.getElementsByTagName('chosen').item(0);
                thisListItem.setAttribute('runName',dataName);
                xmlwrite('Mitodix.config',xmlFile);
            end
            
            % add to matlab paths:
            addpath(genpath('.'));
            
            % Initialize hash-tables:
            self.frames = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
                        
            % Set configuration and logger
            self.Configure = Config();
            self.dataName = self.Configure.DataName;
            self.Log = SingleInstance.Logger;

            % Set debug mode:
            if(~exist('debugMode', 'var') || isempty(debugMode))
                debugMode = 0;
            end
            self.DEBUG = debugMode;
            
            % Creation time:
            self.creationTime = datestr(now,'yyyy-mm-dd_HH-MM-SS-FFF');
            
            % Create DEBUG folder:
            if(self.DEBUG)
                debugFolder = fullfile(self.OutputPath,'DEBUG');
                if(~exist(debugFolder, 'dir'))
                    mkdir(debugFolder);
                end
                
                self.run_log = fullfile(debugFolder, 'MitodixRun.log');
                fid = fopen( self.run_log, 'wt' );
                if(fid == -1)
                    self.run_log = [];
                else
                    fclose(fid);
                end
            else
                self.run_log = [];
            end
        end
        
        % Get Frames list
        function Frames = get.Frames(self)
            Frames = [];
            T = self.FramesNumbers;
            for k=1:length(T)
                t = T(k);
                Frames = [Frames, self.frames(t)];
            end
        end

        % Get cell distances
        function Dists = get.Dists(self)
            if(isempty(self.cells_dists))
                T = self.FramesNumbers;
                for t=T
                    self.cells_dists = [self.cells_dists; self.frames(t).CellsData.NeighborsDists];
                end
            end
            
            Dists = self.cells_dists;
        end
        
        % Get cell drifts
        function Drifts = get.Drifts(self)
            if(isempty(self.cells_drifts))
                T = self.FramesNumbers;
                for kk=2:length(T)
                    [~,drifts] = knnsearch(self.frames(T(kk)).CellsData.Centers,self.frames(T(kk-1)).CellsData.Centers);
                    self.cells_drifts = [self.cells_drifts; drifts];
                end
            end
            
            Drifts = self.cells_drifts;
        end
        
        % Get path of input gray level images
        function InputImPath = get.InputImPath(self)
            InputImPath = self.Configure.InputImPath;
        end
        
        % Get file extention of input gray level images
        function InputImExtension = get.InputImExtension(self)
            InputImExtension = self.Configure.InputImExtension;
        end
        
        % Get file expression of input gray level images
        function InputImExpression = get.InputImExpression(self)
            InputImExpression = self.Configure.InputImExpression;
        end
        
        % Get file expression of input segmentation images
        function InputSegExpression = get.InputSegExpression(self)
            InputSegExpression = self.Configure.InputSegExpression;
        end
        
        % Get path of input segmentation images
        function InputSegPath = get.InputSegPath(self)
            InputSegPath = self.Configure.InputSegPath;
        end
        
        % Get file extention of input segmentation images
        function InputSegExtension = get.InputSegExtension(self)
            InputSegExtension = self.Configure.InputSegExtension;
        end
        
        % Get output folder path for this specific run (according to
        % Mitodix object creation time)
        function OutputPath = get.OutputPath(self)
            if(~isempty(self.Configure.OutputPath))
                % Make sure data set's output folder exists:
                if(~exist(self.Configure.OutputPath, 'dir'))
                    mkdir(self.Configure.OutputPath);
                end
            end
            outputPath = fullfile(self.Configure.OutputPath, ['Results_',self.creationTime]);
            if(~exist(outputPath, 'dir'))
            	mkdir(outputPath);
            end
            OutputPath = outputPath;
        end
        
        % Get data set name
        function DataName = get.DataName(self)
            DataName = self.dataName;
        end

        % Get Ground truth *.mat file path
        function GTPath = get.GTPath(self)
            GTPath = self.Configure.GTPath;
        end
        
        % Get Ground truth
        function GroundTruth = get.GroundTruth(self)
            if(~exist(self.GTPath, 'file'))
                GroundTruth = [];
            else
                GroundTruth = load(self.GTPath);
                GroundTruth = GroundTruth.GT;
                if(self.downscaleFactor ~= 1)
                    GroundTruth.MotherCOM = double(GroundTruth.MotherCOM)/self.downscaleFactor;
                    GroundTruth.Daughter1COM = double(GroundTruth.Daughter1COM)/self.downscaleFactor;
                    GroundTruth.Daughter2COM = double(GroundTruth.Daughter2COM)/self.downscaleFactor;
                end
                
                if(self.dt ~= 1)
                    GroundTruth.BirthFrame = ceil(GroundTruth.BirthFrame/self.dt);
                end
                
                GroundTruth = GroundTruth(GroundTruth.BirthFrame<=self.LastFrameNum, :);
            end
        end
        
        % Get minimal cell length to consider [pixels]
        function MinCellLength = get.MinCellLength(self)
            MinCellLength = self.Configure.MinCellLength;
        end
        
        % Get last frame in Mitodix collection (Frame object)
        function LastFrame = get.LastFrame(self)
            if(self.LastFrameNum)
                LastFrame = self.frames(self.LastFrameNum);
            else
                LastFrame = Frame.empty;
            end
        end
        
        % Get all frame numbers in frames collection
        function FramesNumbers = get.FramesNumbers(self)
        % Get a vector of all frame numbers from frames:
            FramesNumbers = cell2mat(keys(self.frames));
        end
        
        % Get the first frame number in Mitodix collection
        function FirstFrameNum = get.FirstFrameNum(self)
            if(isempty(self.frames))
                FirstFrameNum = 0;
            else
                FirstFrameNum = min(self.FramesNumbers);
            end
        end
        
        % Get last frame number in Mitodix collection
        function LastFrameNum = get.LastFrameNum(self)
            if(isempty(self.frames))
                LastFrameNum = 0;
            else
                LastFrameNum = max(self.FramesNumbers);
            end
        end
       
        % Read frame input images into Mitodix collection
        function Z_dim = ReadFrames(self, startFrame, endFrame, startSeg, rescaleFlag, voxel_idx)
        % Read frame input images into Mitodix collection
        % Signature:
        %   ReadFrames(self, startFrame, endFrame, startSeg)
        % Inputs:
        %   self = Mitodix object.
        %   startFrame = (Default = startSeg+1) The #frame no. to start from.
        %	endFrame = (Default = all frames in input folder) The #frame no. to end at.
        %   startSeg = (Default = 0) the first segmented frame is startSeg+1 (e.g. segmentations are
        %                from frame number 3, startSeg is 2).
        %%
            self.logMsg(sprintf('Read input data...'));
            TIME = tic();
            
            if(~exist('startSeg', 'var') || isempty(startSeg))
                startSeg = 0;
            end

            if(~exist('endFrame', 'var') || isempty(endFrame))
                endFrame = 0;
            end

            if(~exist('startFrame', 'var') || isempty(startFrame))
                startFrame = startSeg+1;
            end
            
            if(~exist('rescaleFlag', 'var') || isempty(rescaleFlag))
                rescaleFlag = ~isequal(self.DataName(1:3), 'DIC');
            end
                
            if(~exist('voxel_idx', 'var') || isempty(voxel_idx))
                voxel_idx = 1;
            end
            
            [grayScaleList, startFrame, endFrame] = self.GetFrameImages(startFrame, endFrame);
            [labeldList, startFrame2, endFrame2] = self.GetFrameImages(startFrame-startSeg, endFrame-startSeg, 1);
            if(~isempty(labeldList))
                endFrame = min(endFrame, endFrame2+startSeg);
                startFrame = max(startFrame, startFrame2+startSeg);
            end
            
            N = endFrame-startFrame+1;
            imPath = self.InputImPath;
            lblPath = self.InputSegPath;
 
%             h = waitbar(0);
            
            fNums = 1:self.dt:N;
            frames_list = cell(length(fNums), 1);
            max_val = -inf;
            min_val = inf;
            Z_dim = 1;
            downscaleFactor = self.downscaleFactor;
            if(downscaleFactor<0)
                self.downscaleFactor = 1;
            end
            minCellLength = round(double(self.Configure.MinCellLength)/self.downscaleFactor);
            parfor n = 1:length(fNums)
                t = n + startFrame - 1;
%                 waitbar(n/N, h, ['Frame ', num2str(t), ' out of ', num2str(endFrame)]);
                if(isempty(labeldList))
                    seg_name = grayScaleList{fNums(n)};
                else
                    seg_name = labeldList{fNums(n)};
                end
                [ image, labeled, Z_dim ] = Frame.Read2D3D( fullfile(imPath, grayScaleList{fNums(n)}), fullfile(lblPath, seg_name), [], rescaleFlag, downscaleFactor );
                
                % Create a Frame object:
                frames_list{n} = Frame(t, image, labeled, grayScaleList{fNums(n)}, [], seg_name, rescaleFlag, minCellLength);
                max_val = max(max_val, double(frames_list{n}.MaxVal));
                min_val = min(min_val, double(frames_list{n}.MinVal)); 
            end
            
            cellsCount = 0;
            for n=1:length(fNums)
                frames_list{n}.RescaleGrayLevel((2^8-1)./(max_val-min_val), rescaleFlag, -min_val);
                cellsCount = cellsCount + length(frames_list{n}.CellsData.Labels);
                self.AddFrame(frames_list{n});
            end
            self.logMsg(sprintf('%d cells, %d frames ', cellsCount, length(fNums)));
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
            
            % Save checkpoints:
            if(self.DEBUG)
                self.SaveCheckpoint();
                framesPath = fullfile(self.OutputPath, 'Frames');
            end
%             delete(h);
            
        end
        
        % Save Mitodix Checkpoint to a file
        function SaveCheckpoint(self)
        % Save Mitodix Checkpoint to a file
        % Signature:
        %   SaveCheckpoint(self, folderPath)
        % Inputs:
        %   self = Mitodix object.
        %%  
            checkpointsPath = fullfile(self.OutputPath, 'Checkpoints');
            if(~exist(checkpointsPath, 'dir'))
                mkdir(checkpointsPath);
            end
            
            mtdx = self;
            timeTag = datestr(now,'yyyy-mm-dd_HHMMSS-FFF');
            try
                save(fullfile(checkpointsPath, ['mtdx_',timeTag,'.mat']), 'mtdx');
            end
         end
         
        % Save Mitodix Frame objects to a file
        function SaveFrames(self, folderPath)
        % Save Mitodix Frame objects to a file
        % Signature:
        %   SaveFrames(self, folderPath)
        % Inputs:
        %   self = Mitodix object.
        %   folderPath = (Default = object's Checkpoints OutputPath). Folder to save
        %                frames *.mat file to.
        %%
            if(~exist('folderPath', 'var') || isempty(folderPath))
                folderPath = fullfile(self.OutputPath, 'Checkpoints');
            end
            if(~exist(folderPath, 'dir'))
                mkdir(folderPath);
            end
            
            frames = self.Frames;
            if(~isempty(frames))
                save(fullfile(folderPath, 'frames.mat'), 'frames');
            end
        end
        
        % Reset Mitodix object and load Frames from a frames.mat file
        function LoadFrames(self, folderPath)
        % Reset Mitodix object and load Frames from a frames.mat file
        % Signature:
        %   LoadFrames(self, folderPath)
        % Inputs:
        %   self = Mitodix object.
        %   folderPath = (Default = object's OutputPath). Folder to read
        %                frames.mat file from.
        %%
            if(~exist('folderPath', 'var') || isempty(folderPath))
                folderPath = self.OutputPath;
            end
            filePath = fullfile(folderPath, 'frames.mat');
            
            if(~exist(filePath, 'file'))
                self.Log.mlog = {LogType.WARNING ...
                                  ,mfilename('class') ...
                                  ,[self.Log.Me,'Cannot load frames from file: ', filePath, ' Does not exist.']};
                return;
            end
            
            self.Reset();
            framesList = load(fullfile(folderPath, 'frames.mat'));
            framesList = framesList.frames;
            N = length(framesList);
            for n=1:N
                self.AddFrame(framesList(n));
            end
        end
         
        % Reset Mitodix object for a new run
        function Reset(self)
        % Reset Mitodix object for a new run
        %%
            % save previous results before reset:
            self.frames = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
            self.creationTime = datestr(now,'yyyy-mm-dd_HHMMSS-FFF');
            
            % Set configuration:
            self.Configure = Config();
            self.dataName = self.Configure.DataName;

            if(~isempty(self.Configure.OutputPath))
                % Make sure data set's output folder exists:
                if(~exist(self.Configure.OutputPath, 'dir'))
                    mkdir(self.Configure.OutputPath);
                end
            end
            
            if(self.DEBUG)
                debugFolder = fullfile(self.OutputPath,'DEBUG');
                if(~exist(debugFolder, 'dir'))
                    mkdir(debugFolder);
                end

                self.run_log = fullfile(mtdx.OutputPath, 'MitodixRun.log');
                fid = fopen( log_file, 'wt' );
                if(fid == -1)
                    self.run_log = [];
                else
                    fclose(fid);
                end
            else
                self.run_log = [];
            end
        end
        
        % Add Frame object to the frames collection
        function AddFrame(self, frame)
        % Add Frame object to the frames collection
        %%
            % Check that input Frame object is not empty:
            if isempty(frame)
                self.Log.mlog = {LogType.WARNING ...
                                  ,mfilename('class') ...
                                  ,[self.Log.Me,' Input Frame is empty.']};
                return;
            end
            
            % If frames already have this frame, notify it is replaced with the new one:
            if(self.frames.isKey(frame.Number))
                self.Log.mlog = {LogType.DEBUG ...
                                  ,mfilename('class') ...
                                  ,[self.Log.Me,' Frame object of frame ', num2str(frame.Number),' is replaced.']};
            end
            
            % Calculate CellDrift:
            prev = frame.Number-1;
            if(self.frames.isKey(prev))
                prevFrame = self.frames(prev);
                frame.SetCellDrift(prevFrame.CellsData.Centers);
            end
            % Add the new frame to the collection:
            self.frames(frame.Number) = frame;
                        
            framesPath = fullfile(self.OutputPath, 'Frames');
            if(~exist(framesPath, 'dir'))
                mkdir(framesPath);
            end
            save(fullfile(framesPath, ['frame', num2str(frame.Number), '.mat']), 'frame');
        end
        
        % Get Frame object from the frames collection
        function frame = GetFrame(self, frameNum)
        % Get Frame object from the frames collection
        %%
        % get Frame of frameNum from frames table:
            frame = Frame.empty;
            if(self.frames.isKey(frameNum))
                frame = self.frames(frameNum);
            end
        end
        
        % Get a list of the frame image files
        function [frameNamesList, startFrame, endFrame] = GetFrameImages(self, startFrame, endFrame, isLabels)
        % Get a list of the frame image files.
        % Inputs:
        %   dataName: data set name in the Mitodix.config.
        %   startFrame: The #frame no. to start from (Defalt=1).
        %   endFrame: The #frame no. to end at (Default: all the images).
        %   isLabels: flag 1=segmentations, 0=intensities (Default=0).
        %%
            errMsg = '';
            
            try
                if(~exist('isLabels', 'var') || isempty(isLabels))
                    isLabels = 0;
                end
                
                % Check if input folder exists:
                if(isLabels)
                    inputPath = self.InputSegPath;
                    ext = self.InputSegExtension;
                    expr = self.InputSegExpression;
                else
                    inputPath = self.InputImPath;
                    ext = self.InputImExtension;
                    expr = self.InputImExpression;
                end
                    
                if(~exist(inputPath, 'dir'))
                    errMsg = ['Frames input folder ', inputPath, ' does not exist.'];
                    error(errMsg);
                end
                
                if isempty(expr)
                    expr = '\w+(\d+).*';
                end

                % Find all image files in input folder:
                fileList = dir(fullfile(inputPath,['*.',ext]));
                fileNames = {fileList.name};
                [t, frameNamesList] = Frame.GetFileNameToken( fileNames, expr );
                N = length(frameNamesList);

                % Check if input folder is not empty:
                if(~N)
                    errMsg = ['There are no relevant files in frames input folder ', inputPath,'.'];
                    warning(errMsg);
                end
                
                if(nargin < 2)
                    startFrame = 1;
                else
                    startFrame = max(startFrame, 1);
                end
                
                if((nargin < 3) || (endFrame < startFrame))
                    endFrame = N;
                else
                    endFrame = min(endFrame, N);
                end
                frameNamesList = frameNamesList(startFrame:endFrame);
                
            catch EX
                if(~isempty(errMsg))
                    self.Log.mlog = {LogType.ERROR ...
                                    ,mfilename('class') ...
                                    ,[self.Log.Me,' Problem getting input frame images files list: ', errMsg]};
                end
                rethrow(EX);
            end
        end
        
        % Write messege into run log
        function logMsg(self, str_msg)
            fprintf(str_msg);
            if(~isempty(self.run_log))
                fid = fopen( self.run_log, 'a' );
                fprintf(fid, str_msg);
                fclose(fid);
            end
        end

        % Create an output video, mark mitosis candidates
        function CreateOutputVideo(self, candidates, pSize, frameRate, startFrame, endFrame)
        % Create an output video, mark mitosis candidates
        % Signature:
        %   CreateOutputVideo(self, candidates, pSize, frameRate, startFrame, endFrame)
        % Inputs:
        %   self: Mitodix object.
        %   candidates: list of candidates for mitosis annotation on video.
        %   pSize: used for margin removal (margin=(pSize-1)/2).
        %   frameRate: Desired frame rate [sec] (Default: 1 sec).
        %   startFrame: initial frame number in video (Default: first frame in
        %   mitodix collection).
        %   endFrame: last frame number in video (Default: last frame in
        %   mitodix collection).
        %%
        
            self.logMsg(sprintf('Video output...'));
            TIME = tic();
            %% Set default input values:
            
            if(~exist('startFrame', 'var') || isempty(startFrame))
                startFrame = self.FirstFrameNum;
            end
            
            if(~exist('endFrame', 'var') || isempty(endFrame))
                endFrame = self.LastFrameNum;
            end
            
            if(~exist('pSize', 'var') || isempty(pSize))
                w = 0;
            else
                w = (pSize-1)/2;
            end
            
            %% Annotate frame by frame:
            vidWriter = VideoWriter.empty;
            mInd = 1;
            dInd = 1;
            
            % Output folder:
            vidFolder = fullfile(self.OutputPath, 'Visualize');
            if(~exist(vidFolder, 'dir'))
                mkdir(vidFolder);
            end

            % choose a diferent color for each candidate:
            colorsNames = {'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'white'};
            colorsRGB = 255*[0, 0, 1; 0, 1, 0; 1, 0, 0; 0, 1, 1; 1, 0, 1; 1, 1, 0; 1, 1, 1];
            colorsN = length(colorsNames);
            N = endFrame - startFrame + 1;
            if(~isempty(candidates))
                candM = candidates(mInd);
                candD = candidates(dInd);
            end
            % Write frames into video:
            for n = 1:N
                t = n + startFrame - 1;
                frame = self.GetFrame(t);
                I = frame.Image;

                if(frame.Z_dim > 1)
                    % TO DO: create 3D video
                    delete(fullfile(vidFolder, ['Results_',self.DataName,'.avi']));
                    return;
                end
                % read input frame image (gray scale):
                
                if(~isa(I, 'uint8'))
                    I = double(I);
                    I = (2^8-1) * I/max(I(:));
                    I = uint8(I);
                end
                I = cat(length(size(I))+1,I,I,I);

                % mothers in this frame:
                while(mInd < length(candidates) && candM.FrameNum < t+1) 
                    mInd = mInd + 1;
                    candM = candidates(mInd);
                end
                while((mInd <= length(candidates)) && (candM.FrameNum == t+1))
                    mSerialNumber = num2str(mInd);

                    colorIdx = mod(mInd, colorsN) + 1;
                    color = colorsRGB(colorIdx, :);
                    colorName = colorsNames(colorIdx);

                    % Mark COM:
                    center = candM.COM;
                    if(frame.Z_dim == 1)
                        zz = 1;
                    else
                        zz = center(end);
                    end
                    
                    if(frame.Z_dim > 1)
                        C = reshape(I(:, :, zz, :),  size(I,1), size(I,2), 3);
                        I(:,:,zz,:) = insertMarker(C, center(1:2) , 'x', ...
                                    'Color', colorName);
                    else
                        I = insertMarker(I, center(1:2) , 'x', ...
                                    'Color', colorName);
                    end
                    
                    % Mark candidate serial number:
                    textLocation = int32([center(1), center(2)-5]);
                    if(frame.Z_dim > 1)
                        C = reshape(I(:, :, zz, :),  size(I,1), size(I,2), 3);
                        I(:,:,zz,:) = insertText(C,textLocation,mSerialNumber,'FontSize',20,'TextColor',colorName,'BoxOpacity',0);
                    else
                        I = insertText(I,textLocation,mSerialNumber,'FontSize',20,'TextColor',colorName,'BoxOpacity',0);
                    end

                    % Mark cell contour:
                    mLabel = candM.MotherLabel;
                    if((~isempty(frame)) && (mLabel~=0))
                        idxs = CalcContour(frame.Labeled == mLabel)~=0;
                        if(frame.Z_dim > 1)
                            IR = I(:,:,:,1); IG = I(:,:,:,2); IB = I(:,:,:,3);
                        else
                            IR = I(:,:,1); IG = I(:,:,2); IB = I(:,:,3);
                        end
                        IR(idxs) = color(1); IG(idxs) = color(2); IB(idxs) = color(3);
                        I = cat(length(size(I)),IR,IG,IB);
                    end
                    mInd = mInd + 1;
                    if(mInd <= length(candidates))
                        candM = candidates(mInd);
                    else
                        candM = [];
                    end
                end

                % daughters in this frame:
                while(dInd < length(candidates) && candD.FrameNum < t) 
                    dInd = dInd + 1;
                    candD = candidates(dInd);
                end
                while((dInd <= length(candidates)) && (candD.FrameNum == t))
                    dSerialNumber = num2str(dInd);

                    colorIdx = mod(dInd, colorsN) + 1;
                    color = colorsRGB(mod(dInd, colorsN) + 1, :);
                    colorName = colorsNames(colorIdx);

                    lDaughter1 = candD.DaughtersLabels(1);
                    lDaughter2 = candD.DaughtersLabels(2);
                    center1 = candD.Daughter1COM;
                    center2 = candD.Daughter2COM;
                    midPoint = round(0.5*(center1 + center2));
                    if(frame.Z_dim == 1)
                        zz = 1;
                    else
                        zz = max(min(round(midPoint(end)), size(I, 3)), 1);
                    end
                    
                    % mark COMs:
                    textLocation = int32(midPoint(1:2));
                    if(frame.Z_dim > 1)
                        C = reshape(I(:, :, zz, :),  size(I,1), size(I,2), 3);
                        I(:,:,zz,:) = insertMarker(C, center1(1:2) , 'x', ...
                                        'Color', colorName);
                        C = reshape(I(:, :, zz, :),  size(I,1), size(I,2), 3);
                        I(:,:,zz,:) = insertMarker(C, center2(1:2) , 'x', ...
                            'Color', colorName);

                        % Mark candidate serial number:
                        C = reshape(I(:, :, zz, :),  size(I,1), size(I,2), 3);
                        I(:,:,zz,:) = insertText(C,textLocation,dSerialNumber,'FontSize',20,'TextColor',colorName, 'BoxOpacity',0);
                    else
                        I = insertMarker(I, center1(1:2) , 'x','Color', colorName);
                        I = insertMarker(I, center2(1:2) , 'x','Color', colorName);
                        % Mark candidate serial number:
                        I = insertText(I,textLocation,dSerialNumber,'FontSize',20,'TextColor',colorName, 'BoxOpacity',0);
                    end
                    
                    % Mark cell contour:
                    if((~isempty(frame)) && (lDaughter1 ~= 0) && (lDaughter2 ~= 0))
                        idxs = (CalcContour(frame.Labeled == lDaughter1) ~= 0) | (CalcContour(frame.Labeled == lDaughter2) ~= 0);
                        if(frame.Z_dim > 1)
                            IR = I(:,:,:,1); IG = I(:,:,:,2); IB = I(:,:,:,3);
                        else
                            IR = I(:,:,1); IG = I(:,:,2); IB = I(:,:,3);
                        end
                        IR(idxs) = color(1); IG(idxs) = color(2); IB(idxs) = color(3);
                        I = cat(length(size(I)),IR,IG,IB);
                    end

                    dInd = dInd + 1;
                    if(dInd <= length(candidates))
                        candD = candidates(dInd);
                    else
                        candD = [];
                    end
                end

                % Add frame into video:
                Iclip = I(w+1:end-w,w+1:end-w, :, :);
                if(frame.Z_dim == 1)
                    imwrite(Iclip, fullfile(vidFolder, ['Vis_', frame.Name]));
                else
                    f = figure;
                    num_x = size(Iclip,2);
                    num_y = size(Iclip,1);
                    Y = (1:num_y);
                    X = (1:num_x);
                    for zz=5:size(I, 3)-5
                        C = reshape(Iclip(:, :, zz, :),  size(Iclip,1), size(Iclip,2), 3);
                        surf(X, Y, 3*zz*ones(size(C,1), size(C,2)), C, 'edgecolor', 'none');
                        hold on;
                    end
                    hold off;
                    axis equal
                    saveas(f, fullfile(vidFolder, ['Vis_', frame.Name]));
                    close(f);
                end
            end
                
            try
                if(~exist('frameRate', 'var') || isempty(frameRate) || (frameRate < 1))
                    frameRate = round(20/self.Configure.FrameRate);
                end
            
                % Create a video writer object:
                vidWriter = VideoWriter(fullfile(vidFolder, ['Results_',self.DataName,'.avi']));
                vidWriter.FrameRate = 10;
                repeat = round(vidWriter.FrameRate/frameRate);
                
                % Open the video writer object:
                open(vidWriter);
                
                expr = ['Vis_', self.InputImExpression];
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
                toc(TIME);
                self.Log.mlog = {LogType.ERROR ...
                                  ,mfilename('class') ...
                                  ,[self.Log.Me,' Problem creating video output for ''',self.DataName,''' data:', EX.message]};
                rethrow(EX);
            end
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
        end

        % Mitosis Detection Framework
        function Detection(self, comperToSvd, startFrame, endFrame)
            if(~exist('comperToSvd', 'var'))
               comperToSvd = self.DEBUG;
            end
           
            if(~exist('startFrame', 'var'))
                startFrame = [];
            end
            if(~exist('endFrame', 'var'))
                endFrame = [];
            end
            
            % Read input frames (Intensity + segmentations):
            Z_dim = self.ReadFrames(startFrame, endFrame);
            
            if(Z_dim == 1)
                %2D images:
            else
                %3D images:
            end

            % Patches:
            [allCandidates, pSize, w, cells_count, features, all_daughter_cell_cands, historyNum, data, weights] = self.GetAllCandidates([],[],comperToSvd);
            featuresVec = features(:,1);
            
            % Outliers detection:
            [otsuThr_mother, otsuThr_daughter] = Mitodix.GetOutliers(featuresVec, all_daughter_cell_cands(:,3), self.DEBUG);
            
            % Outliers cells:
            [motherCandidates, daughter_cell_cands, TPRSim, FPRSim, precisionSim] = self.FindOutliers(allCandidates, featuresVec, all_daughter_cell_cands, otsuThr_mother, otsuThr_daughter);
            
            % Compare to SVD:
            if(self.DEBUG && comperToSvd)
                self.CompareToSVD(allCandidates, historyNum, data, weights, FPRSim, TPRSim, precisionSim);
            end
            
            % Candidate daughters:
            [daughtersCandidates, ~, daughtersPoints, motherPoints] = self.CreateCandidates(motherCandidates, pSize, daughter_cell_cands);

            % Mother-Daughters accossiation:
            [candidatesMatrix, mitoticFeatures, mitoticCandidates, motherPoints, motherCandidates, daughtersPoints, daughtersCandidates] = self.AssociateMotherDaughters(daughtersPoints, motherPoints, daughtersCandidates, motherCandidates, pSize);

            % Mitosis clustering & likelihood:
            [mitosisLikelihood, numClusters] = self.MitosisLikelihood(mitoticFeatures, mitoticCandidates);

            % Decision:
            [finalCandidates, finalMitosisLikelihood] = self.Decision(candidatesMatrix, mitosisLikelihood, mitoticFeatures, daughtersPoints, motherPoints, daughtersCandidates, motherCandidates, [], mitoticCandidates, numClusters);

            % Evaluation:
            if(self.DEBUG)
                self.logMsg(sprintf('Evaluation...'));
                TIME = tic();
               
                % Performance:
                [performance, falsePositiveIdxs, falseNegativeIdxs] = evaluation( self, self.GroundTruth, finalCandidates, -finalMitosisLikelihood, cells_count );
                self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
                
                performanceTbl = struct2table(performance{1});
                writetable(performanceTbl, fullfile(self.OutputPath, 'performanceTbl.txt'));
                save(fullfile(self.OutputPath,'falsesAnalysis.mat'), 'falsePositiveIdxs', 'falseNegativeIdxs');

                struct = performance{1};
                fields = fieldnames(struct);
                for kk = 1:length(fields)
                    self.logMsg(sprintf('    %s: %s\n', fields{kk}, struct.(fields{kk})));
                end
            end

            % Video output:
            self.CreateOutputVideo(finalCandidates); %, pSize);
        end
        
        % Patches Analysis
        function [allCandidates, pSize, w, cells_count, features, daughter_cell_cands,historyNum, data, weights] = GetAllCandidates(self, historyNum, w, comperToSvd)
            self.logMsg(sprintf('Patches extraction...'));
            TIME = tic();

            if(~exist('historyNum', 'var'))
            	historyNum = [];
            end
            if(~exist('w', 'var'))
            	w = [];
            end
            if(~exist('comperToSvd', 'var') || isempty(comperToSvd))
                comperToSvd = 0;
            end
            
            [allCandidates, data, weights, pSize, historyNum, features, daughter_cell_cands] = self.calcPatches(historyNum, w, comperToSvd);
            cells_count = length(allCandidates);
            w = (pSize-1)/2;
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
        end
        
        % Outlier cells
        function  [motherCandidates, daughter_cell_cands, TPRSim, FPRSim, precisionSim] = FindOutliers(self, allCandidates, features, daughter_cell_cands, otsuThr_mother, otsuThr_daughter)
        %%
            if(~exist('otsuThr_mother', 'var') || isempty(otsuThr_mother))
                otsuThr_mother = Mitodix.GetOutliers(features(:,1));
            end
            
            if(~exist('otsuThr_daughter', 'var') || isempty(otsuThr_daughter))
                [~, otsuThr_daughter] = Mitodix.GetOutliers([], daughter_cell_cands(:,3));
            end
            
            % --- Mother outliers: ---
            candIdxSim = find(features(:,1) < otsuThr_mother);
            motherCandidates = allCandidates(candIdxSim);
            
            % --- Daughter outliers: ---
            candIdxSim_daughter = find(daughter_cell_cands(:, 3) < otsuThr_daughter);
            all_daughter_cell_cands = daughter_cell_cands;
            daughter_cell_cands = all_daughter_cell_cands(candIdxSim_daughter, :);
                        
            if(self.DEBUG)
                outputFolder = fullfile(self.OutputPath, 'outliers');
                if(~exist(outputFolder, 'dir'))
                    mkdir(outputFolder);
                end
                save(fullfile(outputFolder, 'outliers.mat'), 'motherCandidates', 'candIdxSim', 'otsuThr_mother', 'allCandidates', 'features', 'daughter_cell_cands', 'all_daughter_cell_cands', 'otsuThr_daughter', '-v7.3');
                
                [ TPs, FPs, TPRSim, FPRSim, precisionSim ] = ROC( self, features(:,1), 'sim', allCandidates, self.GroundTruth, [], otsuThr_mother);
                if(size(features, 2) > 1)
                    ROC( self, features(:, 2), 'sim2', allCandidates, self.GroundTruth );
                    f = figure; scatter(features(FPs, 1), features(FPs, 2), '.b'); hold on;
                    scatter(features(TPs, 1), features(TPs, 2), '*r'); hold off;
                    xlabel('patch sim');
                    ylabel('max overlap');
                    legend({'FP', 'TP'}, 'Location', 'Best');
                    saveas(f, fullfile(outputFolder, 'patches_features.jpg'));
                    saveas(f, fullfile(outputFolder,'patches_features.fig'));
                    close(f);
                end
                
                d_cands = [];
                for kk=1:size(all_daughter_cell_cands, 1)
                    cand = [];
                    cand.BirthFrame = all_daughter_cell_cands(kk, 1);
                    cand.MotherCOM = [-1, -1];
                    cand.Daughter1COM = all_daughter_cell_cands(kk, 4:5);
                    cand.Daughter2COM = [-1, -1];
                    d_cands = [d_cands; cand];
                end
                ROC( self, all_daughter_cell_cands(:, 3), 'daughtersSim', d_cands, self.GroundTruth, [], otsuThr_daughter, [], [], 0);
                
                f = figure;
                [hh, vv] = hist(all_daughter_cell_cands(:, 3), 256);
                bar(vv, hh); hold on; plot(otsuThr_daughter*ones(max(hh), 1), 1:max(hh), '--g'); hold off;
                xlabel('patch sim');
                ylabel('cells count');
                legend({'histogram', 'thr'}, 'Location', 'Best');
                saveas(f, fullfile(outputFolder, 'outlier_daughters.jpg'));
                saveas(f, fullfile(outputFolder,'outlier_daughters.fig'));
                close(f);
            else
                TPRSim = [];
                FPRSim = [];
                precisionSim = [];
            end
        end
        
        function CompareToSVD(self, allCandidates, historyNum, patches, weights, FPRSim, TPRSim, precisionSim)
            self.logMsg(sprintf('SVD based outliers detection (for comparison)...'));
            TIME = tic();
            % Compare to SVD
            
            data = 0*patches;
            % normalization:
            for c=1:size(data, 1)
                mask = (weights(c, :)>0);
                minVal = min(patches(c, mask));
                maxVal = max(patches(c, mask));
                data(c, mask) = (patches(c, mask) - minVal) ./ (maxVal - minVal);
                inds = (data(c, :) < 0);
                data(c, inds) = 0;
            end
            
            [x0, mu, V, S, ~, dim] = self.calcSVD(historyNum, data, 30000, 1, 4);
            mse = self.calcPatchCandidatesSVD(allCandidates, x0, mu, data, weights, V, S, historyNum, dim );
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));

            [~, ~, TPRSvd, FPRSvd, precisionSvd ] = ROC( self, -mse, 'SVD', allCandidates, self.GroundTruth);

            X = cell(2, 1); Y = cell(2, 1); names = cell(2, 1);
            names{1} = 'PCA'; names{2} = 'ours';
            X{1} = FPRSvd; X{2} = FPRSim;
            Y{1} = TPRSvd; Y{2} = TPRSim;
            Z{1} = precisionSvd; Z{2} = precisionSim;
            compareROC( X, Y, Z, names, self.OutputPath );
        end
    
        % Daughters candidates
        function [daughtersCandidates, daughtersSim, daughtersPoints, motherPoints] = CreateCandidates(self, initialCandidates, pSize, daughter_cell_cands)
        % Daughters candidates
        %%
            self.logMsg(sprintf('Daughter candidates (triangultion + symmetry)...'));
            TIME = tic();
            if(mod(pSize, 2) == 0)
                pSize = pSize+1;
            end
            
            w = (pSize - 1)/2;
            lastFrame = self.LastFrame;
            Rd_min = lastFrame.GetDistRadius();
            Rm_min = double(min(self.LastFrame.GetDriftRadius(), ceil(mean(self.Drifts) + 1.5*std(self.Drifts))));
            
            N = length(initialCandidates);
            motherPoints = double(zeros(N, 4));
            daughtersPoints = [];%double(zeros(N, 6));
            framesPath = fullfile(self.OutputPath, 'Frames');
            
%             h = waitbar(0);
            for n=1:N
                motherPoints(n,:) = [initialCandidates(n).FrameNum, initialCandidates(n).COM(1:2), initialCandidates(n).MotherLabel];
                
                frame_file = fullfile(framesPath, ['frame', num2str(initialCandidates(n).FrameNum), '.mat']);
                if(~exist(frame_file, 'file'))
                    continue;
                end
                birthFrame = load(frame_file);
                birthFrame = birthFrame.frame;
                if(isempty(birthFrame))
                    continue;
                end
                
                Rm = birthFrame.GetDriftRadius();
                if(Rm<inf)
                    Rm = max(Rm_min, Rm);
                else
                    Rm = Rm_min;
                end
                Rd = max(Rd_min, birthFrame.GetDistRadius());
                
                % Get all potential daughters:
                motherDists = sqrt(((birthFrame.CellsData.NeighborsMidPoints(:, 1) - initialCandidates(n).COM(1)).^2) + ((birthFrame.CellsData.NeighborsMidPoints(:, 2) - initialCandidates(n).COM(2)).^2));
                daughtersDists = double(floor(birthFrame.CellsData.NeighborsDists'-1));
                ind = find((motherDists <= Rm) & (daughtersDists <= Rd));
                if(isempty(ind))
                    continue;
                end
                pairs = birthFrame.CellsData.Neighbors(ind, :);
                mids = double(birthFrame.CellsData.NeighborsMidPoints(ind, 1:2));
                motherDists = motherDists(ind);
                daughtersDists = birthFrame.CellsData.NeighborsDists(ind);
                
                if(isfield(initialCandidates(n),'PatchDaughterLabel') && ~isempty(initialCandidates(n).PatchDaughterLabel))
                    % only those with one of the daughters detected in prev
                    % stage:
                    ind = pairs(:,1) == initialCandidates(n).PatchDaughterLabel(1) | pairs(:,2) == initialCandidates(n).PatchDaughterLabel(1);
                    if(length(initialCandidates(n).PatchDaughterLabel) == 2)
                        ind = ind | pairs(:,1) == initialCandidates(n).PatchDaughterLabel(2) | pairs(:,2) == initialCandidates(n).PatchDaughterLabel(2);
                    end
                    pairs = pairs(ind, :);
                    mids = mids(ind, :);
                    motherDists = motherDists(ind);
                    daughtersDists = daughtersDists(ind);
                end
                
                for p=1:size(pairs, 1)
                    point = [double(initialCandidates(n).FrameNum), mids(p, 1:2), sort(pairs(p, 1:2)), daughtersDists(p)];
                    if(exist('daughter_cell_cands', 'var') && ~isempty(daughter_cell_cands))
%                         if(N < 500) %low statistics
                            % daughters must be an outlier:
                            ind1 = (daughter_cell_cands(:,1)==point(1)) & (daughter_cell_cands(:,2)==point(4));
                            ind2 = (daughter_cell_cands(:,1)==point(1)) & (daughter_cell_cands(:,2)==point(5));
                            if(~any(ind1) || ~any(ind2))
                                continue;
                            end
%                         else
%                             % One of the daughters must be an outlier:
%                             ind = (daughter_cell_cands(:,1)==point(1)) & (((daughter_cell_cands(:,2)==point(4))) | (daughter_cell_cands(:,2)==point(5)));
%                             if(~any(ind))
%                                 continue;
%                             end
%                         end
                    end
                    daughtersPoints = [daughtersPoints; point];
                end
            end
            daughtersPoints = unique(daughtersPoints, 'rows');
            Nd = size(daughtersPoints, 1);
            daughtersPoints = [daughtersPoints, zeros(Nd, 1)];
            daughtersSim = double(zeros(Nd,1));
            daughtersCandidates = cell(Nd, 1);
            framesPath = fullfile(self.OutputPath, 'Frames');
            parfor n=1:Nd
%                 waitbar(n/Nd, h, ['Candidate ', num2str(n), ' out of ', num2str(Nd)]);
                
                point = daughtersPoints(n, :);
                cand = [];
                
            	cand.FrameNum = point(1);
                
                frame_file = fullfile(framesPath, ['frame', num2str(cand.FrameNum), '.mat']);
                if(~exist(frame_file, 'file'))
                    continue;
                end
                birthFrame = load(frame_file);
                birthFrame = birthFrame.frame;
                if(isempty(birthFrame))
                    continue;
                end
                
                cand.DaughtersMid = point(2:3);
                cand.DaughtersLabels = point(4:5);
                cand.Rd = point(6);
                [ cand.DaughtersIm, cand.DaughtersWeightsIm,  cand.DaughtersMaskIm, cand.DaughterImages, coms, deg ] = birthFrame.GetCellsPatch( cand.DaughtersLabels, pSize, 1, 1 );
                if(isempty(cand.DaughtersIm) || isempty(cand.DaughterImages))
                    continue;
                end
                cand.Daughter1COM = coms(1,:);
                cand.Daughter2COM = coms(2,:);
                cand.orientation = deg;
                cand.DaughtersSim = Frame.PatchSim(cand.DaughterImages{1}.Im, cand.DaughterImages{1}.MaskIm, cand.DaughterImages{2}.Im, cand.DaughterImages{2}.MaskIm);
                daughtersSim(n) = cand.DaughtersSim;
                daughtersCandidates{n} = cand;
                point(7) = cand.DaughtersSim;
                daughtersPoints(n, :) = point;
            end
            rmInd = (cellfun(@isempty,daughtersCandidates));
            daughtersSim(rmInd) = [];
            daughtersPoints(rmInd, :) = [];
            daughtersCandidates = daughtersCandidates(~cellfun(@isempty,daughtersCandidates));
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
            
            if(self.DEBUG)
                outputFolder = fullfile(self.OutputPath, 'Symmetry');
                if(~exist(outputFolder, 'dir'))
                    mkdir(outputFolder);
                end
                save(fullfile(outputFolder, 'Symmetry.mat'), 'daughtersCandidates', 'daughtersSim', 'daughtersPoints', 'motherPoints', 'pSize', '-v7.3');
            end
%             close(h);
        end
        % Mother-daughters association
        function [candidatesMatrix, mitoticFeatures, mitoticCandidates, motherPoints, motherCandidates, daughtersPoints, daughtersCandidates] = AssociateMotherDaughters(self, daughtersPoints, motherPoints, daughtersCandidates, motherCandidates, pSize)
        % Mother-daughters association
        %%
            self.logMsg(sprintf('Candidate Triplets...'));
            TIME = tic();
            if(mod(pSize, 2) == 0)
                pSize = pSize+1;
            end
            w = (pSize - 1)/2;
            
            timeDis = [];
            symm = [];
            distsRatio = [];
            overlapRatio = [];
            candidatesMatrix = [];
            mitoticCandidates = [];
            mitoticFeatures = [];

            Nm = length(motherCandidates);
            Nd = length(daughtersCandidates);
            smallStatFlag = Nm < 10;
            
            if((~Nm) || (~Nd))
                toc(TIME);
                return;
            end
            
            % sort by birthFrame number:
            [motherPoints, idx] = sortrows(motherPoints, 1);
            motherCandidates = motherCandidates(idx);
            [daughtersPoints, idx] = sortrows(daughtersPoints, 1);
            daughtersCandidates = daughtersCandidates(idx);
            
            %symmetric daughters:
            symm_vec = double(daughtersPoints(:,end));
%             MED = median(symm_vec);
%             MAD = median(MED - symm_vec(symm_vec<MED));
%             symmetryThr = max(min(0.8, MED - 3*MAD), 0.4); 
            
            [symm_h, vv] = hist(double(daughtersPoints(:,end)), 0:1/16:1);
            
            if(Nd<50)
                % no statistics
                symmetryThr = 0.1;
            else
                quantized = imquantize(symm_vec, vv(2:end-1), vv(1:end-1));
                symmetryThr = vv(1) + graythresh(quantized)*(vv(end)-vv(1));
                symmetryThr = max(min(0.8, abs(symmetryThr - std(symm_vec))), 0.4);
            end
%             
%             [symm_h, vv] = hist(double(daughtersPoints(:,end)), 256);
%             symmetryThr = vv(1) + graythresh(daughtersPoints(:,end))*(vv(end)-vv(1));
%             symmetryThr = min(0.8, symmetryThr);
            
            symm_cond = (daughtersPoints(:,end) >= symmetryThr);

            Rm_min = double(min(self.LastFrame.GetDriftRadius(), ceil(mean(self.Drifts) + 1.5*std(self.Drifts))));
                
            h = waitbar(0);
            for m=1:Nm
                waitbar(m/Nm, h, ['Candidate ', num2str(m), ' out of ', num2str(Nm)]);
                birthFrameNum = motherPoints(m,1);
                frameDaughtersIdx = find((daughtersPoints(:,1) == birthFrameNum) & symm_cond);
                if(isempty(frameDaughtersIdx))
                    continue;
                end
                birthFrame = self.GetFrame(birthFrameNum);
                dc = mean(double(birthFrame.Image(:)));
                
                Rm = double(birthFrame.GetDriftRadius());
                if(Rm<inf)
                    Rm = max(Rm_min, Rm);
                else
                    Rm = Rm_min;
                end
                
%                 if(smallStatFlag)
%                     Rm = 2*Rm;
%                 end
                
                frameDaughters = daughtersPoints(frameDaughtersIdx,:);
                motherDists = sqrt(double(((frameDaughters(:, 2) - motherPoints(m,2)).^2) + ((frameDaughters(:, 3) - motherPoints(m,3)).^2)));
                frameDaughtersIdx = frameDaughtersIdx(motherDists < Rm);
                motherDistsRatio = double(motherDists(motherDists<Rm));
                
                if(isempty(frameDaughtersIdx))
                    continue;
                end
                
                motherFrame = self.GetFrame(birthFrameNum-1);
                
                cdN = length(frameDaughtersIdx);
                dnCandidates = cell(cdN, 1);
                dnCandidatesMatrix = double(zeros(cdN, 3));
                
                mCands = motherCandidates(m);
                if(length(mCands.COM)==3)
                    [lbld, I] = motherFrame.GetPadded(w, 0, round(mCands.COM(3)));
                else
                    [lbld, I] = motherFrame.GetPadded(w, 0);
                end
                I = double(I);
                I = I - mean(I(:)) + dc;
                I(I<0) = 0;
                
                m_com = round(mCands.COM(1:2)) + [w, w];
                mCands.MotherImages.Im = I(m_com(2)-w:m_com(2)+w, m_com(1)-w:m_com(1)+w);
                mCands.MotherImages.MaskIm = (lbld(m_com(2)-w:m_com(2)+w, m_com(1)-w:m_com(1)+w) == mCands.MotherLabel);
                parfor d=1:cdN
                    daughtersIdx = frameDaughtersIdx(d);
                    daughtersCand = daughtersCandidates{daughtersIdx};
                    if(isempty(daughtersCand))
                        dnCandidates{d} = -1;
                        continue;
                    end
                    
                    % dissimilarity between the mother to each of the
                    % candidate daughters:
                    deg = daughtersCand.orientation;
                    MotherImages = [];
                    MotherImages.Im = RigidTransform(mCands.MotherImages.Im, 0, daughtersCand.orientation, 0);
                    MotherImages.MaskIm = RigidTransform(mCands.MotherImages.MaskIm, 0, daughtersCand.orientation, 1);
                    md1 = Frame.PatchSim(MotherImages.Im, MotherImages.MaskIm, daughtersCand.DaughterImages{1}.Im, daughtersCand.DaughterImages{1}.MaskIm, 1);
                    md2 = Frame.PatchSim(MotherImages.Im, MotherImages.MaskIm, flipud(daughtersCand.DaughterImages{2}.Im), flipud(daughtersCand.DaughterImages{2}.MaskIm), 1);
                    candTimeDis = 1 - max(md1, md2);

                    dnCandidatesMatrix(d, :) = [birthFrameNum, m, daughtersIdx];
                    
                    cand = mCands;
                    cand.MotherImages = MotherImages;
                    cand.DaughtersMid = daughtersCand.DaughtersMid;
                    cand.DaughtersLabels = daughtersCand.DaughtersLabels;
                    cand.Rd = daughtersCand.Rd;
                    cand.DaughtersIm = daughtersCand.DaughtersIm;
                    cand.DaughtersWeightsIm = daughtersCand.DaughtersWeightsIm;
                    cand.DaughtersMaskIm = daughtersCand.DaughtersMaskIm;
                    cand.DaughterImages = daughtersCand.DaughterImages;
                    cand.Daughter1COM = daughtersCand.Daughter1COM;
                    cand.Daughter2COM = daughtersCand.Daughter2COM;
                    cand.DaughtersOrientation = daughtersCand.orientation;
                    cand.DaughtersSim = daughtersCand.DaughtersSim;
                    cand.TimeDiss = candTimeDis;
                    cand.MotherDistsRatio = motherDistsRatio(d);
                    cand.OverlapRatio = Frame.PatchSim(cand.MotherImages.Im, cand.MotherImages.MaskIm, cand.DaughtersIm, cand.DaughtersMaskIm);
                    dnCandidates{d} = cand;
                end
                
                for d=1:cdN
                    cand = dnCandidates{d};
                    if(~isstruct(cand) && (cand == -1))
                        continue;
                    end
                    mitoticCandidates = [mitoticCandidates; cand];
                    timeDis = [timeDis; cand.TimeDiss];
                    symm = [symm; cand.DaughtersSim];
                    distsRatio = [distsRatio; cand.MotherDistsRatio];
                    overlapRatio = [overlapRatio; cand.OverlapRatio];
                    candidatesMatrix = [candidatesMatrix; dnCandidatesMatrix(d, :)];
                end
            end
            
            if(smallStatFlag)
                distsRatio = (1-distsRatio./max(distsRatio));
            else
                distsRatio = 1-(distsRatio./Rm);
            end
            motherDaughtersSim = (1-timeDis).*overlapRatio;
            mitoticFeatures = [motherDaughtersSim, symm, distsRatio, timeDis, overlapRatio];
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
            
            if(self.DEBUG)
                outputFolder = fullfile(self.OutputPath, 'Symmetry');
                if(~exist(outputFolder, 'dir'))
                    mkdir(outputFolder);
                end
                
                save(fullfile(outputFolder, 'likelihood.mat'), 'candidatesMatrix', 'mitoticFeatures', 'motherDaughtersSim', 'symm', 'mitoticCandidates', 'pSize', 'symmetryThr', '-v7.3');
                f = figure; bar(vv, symm_h); hold on; plot(symmetryThr*ones(ceil(max(symm_h)), 1), 1:ceil(max(symm_h)), '--g'); title('daughters similarity thr');
                saveas(f, fullfile(outputFolder, 'daughters_similarity_hist.jpg'));
                saveas(f, fullfile(outputFolder, 'daughters_similarity_hist.fig'));
                close(f);
            end
            close(h);
            
        end

        % Mitosis Likelihood
        function [mitosisLikelihood, numClusters] = MitosisLikelihood(self, mitoticFeatures, mitoticCandidates)
            self.logMsg(sprintf('Candidates likelihood (clustering)...'));
            TIME = tic();
            [mitosisLikelihood, gm, c] = softClustering( mitoticFeatures );
            numClusters = length(gm.ComponentProportion);
            features_names = {'mother-daughters similarity', 'daughters symmetry'};
%             if(gm.ComponentProportion(c) > 0.5)
%                 % alternative feature space separation:
%                 mitoticFeatures2 = [mitoticFeatures(:,3), mitoticFeatures(:,2)];
%                 [mitosisLikelihood2, gm2, c2] = softClustering(mitoticFeatures2);
%                 if(gm2.ComponentProportion(c2)<gm.ComponentProportion(c))
%                     mitosisLikelihood = mitosisLikelihood2;
%                     gm = gm2;
%                     c = c2;
%                     mitoticFeatures = mitoticFeatures2;
%                     features_names = {'mother-daughters proximity', 'daughters similarity'};
%                 end
%             end
            
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
            
            if(self.DEBUG && ~isempty(mitoticFeatures))
                mitoticFeaturesSpaceAnalysis( mitoticFeatures, gm, c, self.OutputPath, self.GroundTruth, mitoticCandidates, 5, 0, features_names );
            end
            
        end
        
        % Mitotic Decision (Int. Prog.)
        function [finalCandidates, finalMitosisLikelihood] = Decision(self, candidatesMatrix, mitosisLikelihood, mitoticFeatures, daughtersPoints, motherPoints, daughtersCandidates, motherCandidates, Nm, mitoticCandidates, numClusters)
        % Mitotic Decision
        %%
            orgMitosisLikelihood = mitosisLikelihood;
            if(~exist('Nm', 'var'))
                Nm = -1;
            end
            
            if(~exist('numClusters', 'var') || isempty(numClusters))
                numClusters = 2;
            end
            
            self.logMsg(sprintf('Mitotic Decision (integer programming)...'));
            TIME = tic();
            if(isempty(candidatesMatrix)|| isempty(mitosisLikelihood))
                finalCandidates = []; finalMitosisLikelihood = [];
                toc(TIME);
                return;
            end
            
            fullTimeDiss = mitoticFeatures(:,1);
            
            % Consecutive frames optimization:
            removeCandIdx = self.timeOptimization(candidatesMatrix, mitosisLikelihood, daughtersPoints, motherPoints);
            mitosisLikelihood(removeCandIdx) = -1;
%           
            [hh, vv] = hist(mitosisLikelihood(mitosisLikelihood>0), 256);
            likelihoodThr = max(vv(2), 0.005);
            if(numClusters==1)%(length(mitosisLikelihood) < 50)
                % No statistics:
                [maxH, maxV] = max(hh);
                if((maxV ~= 1) || ((maxH/sum(hh))<0.1))
                    likelihoodThr = 0;
                end
%                 likelihoodThr = vv(1) + graythresh(mitosisLikelihood(mitosisLikelihood>0))*(vv(end)-vv(1));
%                 likelihoodThr = min(0.2, likelihoodThr);
            end
            
            mitosisLikelihood(mitosisLikelihood < likelihoodThr) = -2;
            
            % a priori mitosis probility:
            if(Nm<0)
                Pm = (1/self.Configure.CellCycle);
            end
            
            fNums = unique(candidatesMatrix(:,1));
            fNums = fNums(fNums>0);
            removeCandIdx = 0;
            while(~isempty(removeCandIdx))
                choosenIdx = [];
                for ii=1:length(fNums)
                    t = fNums(ii);
                    tCandIdx = find((candidatesMatrix(:,1) == t) & (mitosisLikelihood >= 0));

                    if(isempty(tCandIdx))
                        continue;
                    end

                    if(length(tCandIdx) == 1)
                        choosenIdx = [choosenIdx; tCandIdx];
                        continue;
                    end

                    mothersIdx = candidatesMatrix(tCandIdx, 2);
                    daughtersIdx = candidatesMatrix(tCandIdx, 3);

                    motherLabels = motherPoints(mothersIdx, 2);
                    daughter1Labels = daughtersPoints(daughtersIdx, 4);
                    daughter2Labels = daughtersPoints(daughtersIdx, 5);
                    mitosisScores = mitosisLikelihood(tCandIdx);
                    tCandidates = table(motherLabels, daughter1Labels, daughter2Labels, mitosisScores);
                    tCandidates.Properties.VariableNames = {'MotherLabel', 'Daughter1Label', 'Daughter2Label', 'MitosisScore'};
                    tCandidates = table2struct(tCandidates);

                    % Probable number of mitosis events in this frame:
                    if(Nm<0)
                        birthFrame = self.GetFrame(t);
                        cellsNum = birthFrame.NumberOfCells;
                        if(cellsNum > 50)
                            % use handreds quantization:
                            cellsNum = 100*ceil(birthFrame.NumberOfCells/100);
                        end
                        Nm = ceil(cellsNum*Pm);
                    end

                    % Solve integer programming problem:
                    solver = IntProg(tCandidates, Nm);
                    [~, idxs] = solver.solve();
                    choosenIdx = [choosenIdx; tCandIdx(idxs)];
                end
                
                finalFeatures = mitoticFeatures(choosenIdx, :);
                finalCandidatesMatrix = candidatesMatrix(choosenIdx, :);
                finalMitosisLikelihood = mitosisLikelihood(choosenIdx);
                timeDiss = fullTimeDiss(choosenIdx);
                finalCandidates = motherCandidates(candidatesMatrix(choosenIdx, 2));
                finalDaughtersCandidates = daughtersCandidates(candidatesMatrix(choosenIdx, 3));
                for n=1:length(choosenIdx)
                    cnd = finalDaughtersCandidates{n};
                    finalCandidates(n).DaughtersMid = cnd.DaughtersMid;
                    finalCandidates(n).DaughtersLabels = cnd.DaughtersLabels;
                    finalCandidates(n).Rd = cnd.Rd;
                    finalCandidates(n).DaughtersIm = cnd.DaughtersIm;
                    finalCandidates(n).DaughtersWeightsIm = cnd.DaughtersWeightsIm;
                    finalCandidates(n).DaughtersMaskIm = cnd.DaughtersMaskIm;
                    finalCandidates(n).DaughterImages = cnd.DaughterImages;
                    finalCandidates(n).Daughter1COM = cnd.Daughter1COM;
                    finalCandidates(n).Daughter2COM = cnd.Daughter2COM;
%                     finalCandidates(n).DaughtersOrientation = cnd.DaughtersOrientation;
                    finalCandidates(n).DaughtersSim = cnd.DaughtersSim;
%                     finalCandidates(n).MotherImages = cnd.MotherImages;
                    finalCandidates(n).TimeDiss = timeDiss(n);
                    finalCandidates(n).MotherDistsRatio = mitoticFeatures(n,3);
                    finalCandidates(n).OverlapRatio = mitoticFeatures(n,4);
                    finalCandidates(n).Likelihood = finalMitosisLikelihood(n);
                end
                
%                 % Consecutive frames optimization:
%                 tmpLikelihood = -3*ones(size(mitosisLikelihood));
%                 tmpLikelihood(choosenIdx) = mitosisLikelihood(choosenIdx);
%                 removeCandIdx = self.timeOptimization(candidatesMatrix, tmpLikelihood, daughtersPoints, motherPoints);
%                 mitosisLikelihood(removeCandIdx) = -1;
                  removeCandIdx = [];
            end
            self.logMsg(sprintf('(%.3f sec)\n', toc(TIME)));
            
            if(self.DEBUG)
                outputFolder = fullfile(self.OutputPath, 'Symmetry');
                if(~exist(outputFolder, 'dir'))
                    mkdir(outputFolder);
                end
                mitosisLikelihood = orgMitosisLikelihood;
                save(fullfile(outputFolder, 'decision.mat'), 'finalCandidates', 'finalMitosisLikelihood', 'mitosisLikelihood', 'likelihoodThr', 'numClusters', '-v7.3');
                
                f = figure; scatter(mitoticFeatures(:, 1), mitoticFeatures(:, 2), '.'); hold on; scatter(finalFeatures(:, 1), finalFeatures(:,2), 'og'); hold off;
                saveas(f, fullfile(outputFolder, 'finalScatter.jpg'));
                saveas(f, fullfile(outputFolder, 'finalScatter.fig'));
                close(f);
                f = figure; scatter3(mitoticFeatures(:, 1), mitoticFeatures(:, 2), mitoticFeatures(:, 3), '.'); hold on; scatter3(finalFeatures(:, 1), finalFeatures(:,2), finalFeatures(:, 3), 'og'); hold off;
                saveas(f, fullfile(outputFolder, 'finalScatter3.jpg'));
                saveas(f, fullfile(outputFolder, 'finalScatter3.fig'));
                close(f);
                try
                    ROC( self, -mitosisLikelihood(mitosisLikelihood>0), 'likelihood', mitoticCandidates(mitosisLikelihood>0), self.GroundTruth, [], -likelihoodThr );
                end
            end
  
        end
        
        % Discard margin cells
        function filtered_GT = FilteredGT(self, w)
        % Discard margin cells.
        %%
            if(exist('w', 'var') && ~isempty(w))
                w = max(w, 50);
            end
            if(~exist('w', 'var') || isempty(w) || (w==0))
                filtered_GT = self.GroundTruth;
            else
                imSizeYX = size(self.LastFrame.Image);
                filtered_GT = RemoveBorderCells( self.GroundTruth, w, imSizeYX );
            end
        end
        
%         function fastCandidtes(self, dt)
%             if(~exist('dt', 'var') || isempty(dt))
%                 dt = 1;
%             end
%             framesPath = fullfile(self.OutputPath, 'Frames');
%             fNums = self.FirstFrameNum:dt:self.LastFrameNum;
%             parfor n=2:length(fNums)
%                 frame_file = fullfile(framesPath, ['frame', num2str(frameNum-1), '.mat']);
%                 frame = load(frame_file);
%                 prevFrame = frame.frame;
%                 
%                 frame_file = fullfile(framesPath, ['frame', num2str(frameNum), '.mat']);
%                 frame = load(frame_file);
%                 curFrame = frame.frame;
%                 
%                 overlaps = prevFrame.Labeled;
%                 overlaps(curFrame.Labeled == 0) = 0;
%             end
%         end
    function [removeCandIdx, removeMotherCandIdx, removeDaughtersCandIdx] = TimeOptimization(self, candidatesMatrix, mitosisLikelihood, daughtersPoints, motherPoints)
        [removeCandIdx, removeMotherCandIdx, removeDaughtersCandIdx] = self.timeOptimization(candidatesMatrix, mitosisLikelihood, daughtersPoints, motherPoints);
    end
    end

    
%% Private methods
    methods (Access = private)
        
        % Patch extraction
        function [allCandidates, data, weights, pSize, historyNum, features, daughter_cell_cands] = calcPatches(self, historyNum, w, comperToSvd)
        % Patch extraction
        %%
            if(~exist('historyNum', 'var') || isempty(historyNum))
            	historyNum = 2;
            end
            if(~exist('comperToSvd', 'var') || isempty(comperToSvd))
            	comperToSvd = 0;
            end
            
            startFrame = self.FirstFrameNum;
            endFrame = self.LastFrameNum;
            
            framesNumbers = self.FramesNumbers;
            %framesNumbers = framesNumbers((framesNumbers >= startFrame) & (framesNumbers > historyNum-1));
            N = length(framesNumbers);
            
            lastFrame = self.GetFrame(endFrame);
            w_m = round(0.5*lastFrame.DistBetweenNeigbors.Mean + 0.5*lastFrame.CellLength.Mean + (historyNum-1)*lastFrame.CellDrift.Mean);
            
            if(~exist('w', 'var') || isempty(w))
                w = 0;
            end
            w = max(w, w_m);
            assert(w>0, 'Patch size must be positive!');
            
            pSize = 2*w+1;
            
            patches_data = cell(N, 1);
            daughter_patches_data = cell(N, 1);
            framesPath = fullfile(self.OutputPath, 'Frames');
            seqEnd =  self.frameNum;
            parfor n=1:N
                [mother_patches_data{n}, daughter_patches_data{n}] = Mitodix.calcInitialPatchFrameCandidates(framesPath, framesNumbers(n), w);
                if(~isnan(seqEnd) && (mod(framesNumbers(n), seqEnd)==0))
                    daughter_patches_data{n} = [];
                end
                if(~isnan(seqEnd) && (mod(framesNumbers(n), seqEnd)==1))
                    mother_patches_data{n} = [];
                end
            end
            
            allCandidates = [];
%             kk = 0;
            data = [];
%             masks = [];
            weights = [];
%             daughterLabeles = [];
%             motherLabels = [];
            features = [];
%             cellsCount = [];
%             patches = [];
            daughter_cell_cands = [];
            h = waitbar(0);
            for n=1:N
                t = framesNumbers(n);
                waitbar(n/N, h, ['Frame ', num2str(t), ' out of ', num2str(endFrame)]);
              
                % --- Mother candidates: ---
                frame_patches = mother_patches_data{n};
                M = length(frame_patches);
                for mm=1:M
                    cand = [];
                    cell_data = frame_patches(mm);
    %                 data = [data; frame_patches.data];
    %                 masks = [masks; frame_patches.masks];
                    
    %                     motherLabels = [motherLabels; frame_patches.motherLabels];
    %                     daughterLabeles = [daughterLabeles; frame_patches.daughterLabeles];
    %                 cellsCount = [cellsCount; frame_patches.cellsCount];
                    features = [features; cell_data.features];
                    
                    if(comperToSvd)
                        data = [data; cell_data.patches];
                        weights = [weights; cell_data.weights];
                    end
                    
                    cand.FrameNum = cell_data.birth_frame;
                    cand.COM = cell_data.centers;%(mm,:);
                    cand.MotherLabel = cell_data.motherLabels;%(mm);
                    cand.PatchSim = cell_data.features;%(mm, :);
                    %cand.MotherImages.Im = cell_data.motherImages;%reshape(frame_patches.motherImages(mm, :), [pSize, pSize]);
                    %cand.MotherImages.MaskIm = cell_data.motherMasks;%reshape(frame_patches.motherMasks(mm, :), [pSize, pSize]);
                    allCandidates = [allCandidates; cand];
                end
                
                % --- Daughter candidates: ---
                frame_patches = daughter_patches_data{n};
                M = length(frame_patches);
                for mm=1:M
                    cell_data = frame_patches(mm);
                    daughter_cell_cands = [daughter_cell_cands; [double(cell_data.birth_frame), double(cell_data.daughterLabels), cell_data.features(:), double(cell_data.centers)]];
                end
            end
            close(h);
            outputFolder = fullfile(self.OutputPath, 'patches');
            if(~exist(outputFolder, 'dir'))
                mkdir(outputFolder);
            end
            if(self.DEBUG)
                save(fullfile(outputFolder, 'patches.mat'), 'allCandidates', 'features', 'data', 'weights', 'pSize', 'historyNum', 'daughter_cell_cands', '-v7.3');%'patches', 'masks', 'cellsCount', '-v7.3');
            end
        end
        
        % Consecutive frames optimization
        function [removeCandIdx, removeMotherCandIdx, removeDaughtersCandIdx] = timeOptimization(~, candidatesMatrix, mitosisLikelihood, daughtersPoints, motherPoints)
        % Consecutive frames optimization
        %%
            removeCandIdx = []; removeMotherCandIdx = []; removeDaughtersCandIdx = [];
            Nm = size(motherPoints, 1);
            for nm=1:Nm
               motherIdx = nm;
               % daughter in frame t-1 is mother of division in frame t: 
               cond1 = (daughtersPoints(:,1) == motherPoints(nm, 1) - 1);
               cond2 = ((daughtersPoints(:,4) == motherPoints(nm, 4)) | (daughtersPoints(:,5) == motherPoints(nm, 4)));               
               daughtersIdx = find(cond1 & cond2);
               if(isempty(daughtersIdx))
                   continue;
               end
               motherCandsIdx = find((candidatesMatrix(:,2) == nm)  & (mitosisLikelihood>0));
               if(isempty(motherCandsIdx))
                   continue;
               end
               daughtersCandsIdx = [];
               for kk=1:length(daughtersIdx)
                   daughtersCandsIdx = [daughtersCandsIdx; find(candidatesMatrix(:,3) == daughtersIdx(kk))];
               end
               likelihoodM = mitosisLikelihood(motherCandsIdx);
               likelihoodD = mitosisLikelihood(daughtersCandsIdx);
               if(max(likelihoodD) >= max(likelihoodM))
                   removeCandIdx = [removeCandIdx; motherCandsIdx];
                   removeMotherCandIdx = [removeMotherCandIdx; nm];
                   motherPoints(motherCandsIdx, 4) = -2;
               else
                   removeCandIdx = [removeCandIdx; daughtersCandsIdx];
                   removeDaughtersCandIdx = [removeDaughtersCandIdx; daughtersIdx];
                   daughtersPoints(daughtersIdx,4) = -1;
                   daughtersPoints(daughtersIdx,5) = -1;
               end
            end
        end
        
        % SVD
        function [x0, mu, V, S, x, dim] = calcSVD(self, historyNum, data, nSamples, usePermutations, scaleFactor)
        % SVD
        %%
            if(~exist('usePermutations', 'var') || isempty(usePermutations))
                usePermutations = 0;
                permutations = [];
            else
                L = 4;
                permutations = cell(L, 1);
            end
            
            nCand = size(data, 1);
            if(~exist('nSamples', 'var') || isempty(nSamples))
                nSamples = nCand;
            end
            
            % Choose subset for PC calculation:
            if(nCand > nSamples)
                rng(1);
                ind = randperm(nCand, nSamples);
            else
                ind = 1:nCand;
            end
            data = data(ind, :);
            nCand = size(data, 1);
            
            sz = size(data, 2) / historyNum;
            pSize = sqrt(sz);
            
            if(~exist('scaleFactor', 'var') || isempty(scaleFactor))
                scaleFactor = 1;
                sampleSize_p = pSize;
            else
                w = (pSize-1)/2;
                sampleSize_p = 2*round(w/scaleFactor)+1;
            end
            
            % permutations:
            if(usePermutations)
                sz_p = sampleSize_p^2; 
                for jj = 1:L
                    permutations{jj} = zeros(nCand, historyNum*sampleSize_p^2);
                end
                
                for c=1:nCand
                    im = reshape(data(c, :), [pSize, pSize*historyNum]);
                    
                    % downsample:
                    if(scaleFactor ~= 1)
                        im = imresize(im, [sampleSize_p, sampleSize_p*historyNum], 'bicubic');
                    end
                    
                    % permutations:
                    permutations{1}(c,:) = im(:);                    
                    for hh=1:historyNum
                        subIm = im(:, sampleSize_p*(hh-1)+1:sampleSize_p*hh);
                        imInd = sz_p*(hh-1)+1:sz_p*hh;
                        im_p = rot90(subIm, 2); permutations{2}(c,imInd) = im_p(:);
                        im_p = flipud(subIm); permutations{3}(c,imInd) = im_p(:);
                        im_p = fliplr(subIm); permutations{4}(c,imInd) = im_p(:);
                    end
                    
                    for jj=1:L
                        % Normalization:
                        imVec = permutations{jj}(c,:);
                        imVec(imVec<0) = 0;
                        imVec = imVec / max(imVec);
                        permutations{jj}(c,:) = imVec;
                    end
                end
                
                x = [];
                for jj = 1:L
                    x = [x; permutations{jj}];
                end
            else
                x = data;
            end
            
            sz = size(x, 2) / historyNum;
            pSize = sqrt(sz);
            dim = historyNum;
            
            x = x/max(x(:));
            
            muVec = mean(x,1);
            mu = repmat(muVec, [size(x, 1), 1]);
            x0 = x -  mu;
            [~, S, V] = svd(x0, 'econ');
            
            if(self.DEBUG)
                outputFolder = fullfile(self.OutputPath, 'SVD');
                if(~exist(outputFolder, 'dir'))
                    mkdir(outputFolder);
                end
                save(fullfile(outputFolder, 'SVD.mat'), 'x0', 'mu', 'V', 'S', 'x', 'dim', '-v7.3');
            end
        end
        
        % SVD-based outliers detection
        function mse = calcPatchCandidatesSVD(self, candidates, x0, mu, data, weights, V, S, historyNum, dim, plotFlag)
        % SVD-based outliers detection
        %%
        
            if(~exist('plotFlag', 'var') || isempty(plotFlag))
                plotFlag = 0;
            end
            
            outputFolder = fullfile(self.OutputPath, 'SVD');
            if(~exist(outputFolder, 'dir'))
                mkdir(outputFolder);
            end
            outputPlots = fullfile(outputFolder, 'plots');
            if(~exist(outputPlots, 'dir'))
                mkdir(outputPlots);
            end

            % Choose number of PCs:
            nCand = length(candidates);
            l = diag(S);
            nobs = size(x0, 1);
            latent = (l.^2)/(nobs-1);
            energy = cumsum(latent)./sum(latent);
            
            latent_s = smooth(smooth(smooth(latent)));
            latent_s = latent - min(latent);
            latent_s = latent_s / max(latent);
            
            d_latent = smooth(smooth(smooth((diff(latent_s)))));
            
            settlingFactor = 0.02;
            latent_err = abs(d_latent - d_latent(end));
            latent_err_tol = settlingFactor*max(latent_err);
            p = find(latent_err > latent_err_tol, 1, 'last') + 1;

            f = figure;
            subplot(2,1,1);
            plot(energy, 'linewidth', 1.5); hold on;
            title('energy'); xlabel('#eigen'); ylabel('energy');
            plot(p, energy(p), 'r*', 'linewidth', 1.5); hold off;
            ax1 = gca;
            
            subplot(2,1,2);
            plot(2:length(d_latent)+1, d_latent, 'linewidth', 1.5); hold on; xlabel('#eigen');
            plot(p, d_latent(p-1), 'r*', 'linewidth', 1.5); hold off;
            legend('latent diff', 'auto calculated p', 'Location', 'south');
            linkaxes([ax1, gca], 'x');
            
            saveas(f, fullfile(outputPlots, 'energy.jpg'));
            saveas(f, fullfile(outputPlots,'energy.fig'));
            
            xlim([1, 50]);
            saveas(f, fullfile(outputPlots, 'energy_zoom.jpg'));
            saveas(f, fullfile(outputPlots,'energy_zoom.fig'));
            
            close(f);

            PC = V(:,1:p);
            l = diag(S);
            S_hat = diag([l(1:p); zeros(length(l)-p,1)]);
            ef = S_hat*V' + mu(1:size(S_hat, 1),:); % eigen cells

            sz_p = (size(x0, 2)/dim);
            ps = sqrt(sz_p);
            
            % plot eigen cells:
            f = figure;
            b = ceil(sqrt(p));
            q = ceil(p/b);
            
            for pp=1:p
                a = ef(pp,:);
                ef_cell = reshape(ef(pp,:), [ps, ps*dim]);
                subplot(b,q,pp); title(['eigen cell #',num2str(p)]); imagesc(ef_cell);
                if(pp == 1)
                    colormap gray;
                    cc = get(gca, 'clim');
                else
                    set(gca, 'clim', cc);
                end
                set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
            end
            saveas(f, fullfile(outputPlots, 'eigen_cells.jpg'));
            saveas(f, fullfile(outputPlots,'eigen_cells.fig'));
            close(f);
            
            % Reconstruction is not measured on permutations:
             sz = size(data, 2)/historyNum;
             pSize = sqrt(sz);
             
             if(sz_p ~= sz)
                 y = zeros(nCand, sz_p*historyNum);
                 for ii=1:nCand
                     downScale = imresize(reshape(data(ii, :), [pSize, pSize*historyNum]), [ps, ps*historyNum], 'bicubic');
                     y(ii, :) = downScale(:);
                 end
             else
                 y = data(1:nCand, :);
             end
             
             if(dim > historyNum)
                 y = [y, y(:,end-pSize^2+1:end) - y(:,1:pSize^2)];
             end

             y = y/max(y(:));
             muY = repmat(mean(y,1), [nCand, 1]);
             y0 = y - muY;

            % Reconstruction error:
             y_rec = y0*(PC*PC') + muY; % reconstruction
             
            mse = zeros(nCand, 1);
            data = data/max(data(:));
            for ii=1:nCand
                rec = y_rec(ii, :);
                org = y(ii, :);
                
                if(sz_p ~= sz)
                    wghts = imresize(reshape(weights(ii, :), [pSize, pSize]), [ps, ps], 'bicubic');
                    wghts = repmat(wghts(:)', [1, dim])';
                else
                    wghts = repmat(weights(ii, :), [1, dim]);
                end
                err = (org - rec)'/max(org);
                mse(ii) = (err'*(wghts.*err))/sum(wghts);
                
                if(self.DEBUG && plotFlag)
                    pSize_p = sqrt(sz_p);
                    
                    f = figure;
                    subplot(4,1,1); imagesc(reshape(data(ii, :), [pSize, pSize*dim])); cc = get(gca, 'clim'); ax1 = gca;
                    set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]); colorbar; title('data');
                    subplot(4,1,2); imagesc(reshape(rec, [pSize_p, pSize_p*dim])); set(gca, 'clim', cc); ax2 = gca;
                    set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]); colorbar; title('reconstruction');
                    subplot(4,1,3); imagesc(reshape(err, [pSize_p, pSize_p*dim])); ax3 = gca;
                    set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]); colorbar; title('error');
                    subplot(4,1,4); imagesc(reshape(wghts, [pSize_p, pSize_p*dim]));
                    set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]); colorbar; title('weights');
                    linkaxes([ax2, ax3, gca]); colormap gray; 
                    close(f);
                    save(fullfile(outputFolder,'SVD_mse.mat'), 'mse', '-v7.3');
                end
            end
        end     
    end
    
    methods (Static)
        % Outliers cells
        function [otsuThr_mother, otsuThr_daughter] = GetOutliers(motherFeaturesVec, daughterFeaturesVec, debugFlag)
            if(~exist('debugFlag', 'var') || isempty(debugFlag))
                debugFlag = 1;
            end
            
            % --- Mother Outliers: ---
            if(exist('motherFeaturesVec', 'var') && ~isempty(motherFeaturesVec))
                otsuThr = AutoThr( motherFeaturesVec, 'sim', debugFlag);
                g = sum(motherFeaturesVec < otsuThr)/length(motherFeaturesVec);
                if(g > 0.1)
                    otsuThr2 = AutoThr( motherFeaturesVec(motherFeaturesVec < otsuThr), 'timeSim', debugFlag);
                    otsuThr = min(otsuThr, otsuThr2);
                end
                otsuThr_mother = otsuThr;
            else
                otsuThr_mother = [];
            end
            
            % --- Daughter outliers: ---
            if(exist('daughterFeaturesVec', 'var') && ~isempty(daughterFeaturesVec))
                otsuThr_daughter = AutoThr( daughterFeaturesVec, 'simDaughter', debugFlag);
                
                g = sum(daughterFeaturesVec < otsuThr_daughter)/length(daughterFeaturesVec);
                if(g > 0.2)
                    otsuThr2 = AutoThr( daughterFeaturesVec(daughterFeaturesVec < otsuThr_daughter), 'simDaughter', debugFlag);
                    otsuThr_daughter = min(otsuThr_daughter, otsuThr2);
                end                
            else
                otsuThr_daughter = [];
            end
        end
        
        function [patch, softMask, mask, centerLabel, status, deg, dc_out] = calcFrameCellPatches(framesPath, frameNum, center, w, dim, refMasked, registerFlag, ref_deg, dc)
            if(~exist('dc', 'var'))
                dc = [];
            end
            frame_file = fullfile(framesPath, ['frame', num2str(frameNum), '.mat']);
            if(~exist(frame_file, 'file'))
                patch = []; softMask = []; mask = [];
                centerLabel = -1;
                status = 0;
                deg = 0;
                dc_out = 0;
                return;
            end
            
            status = 1;
            se = strel('disk',1);
            [optimizer, metric] = imregconfig('monomodal');
            optimizer.MaximumIterations = 1000;
            optimizer.MaximumStepLength = 10^-3;
            
            if(~exist('refMasked', 'var'))
                refMasked = [];
            end

            if(~exist('ref_deg', 'var'))
                ref_deg = [];
            end
            deg = ref_deg;
            
            frame = load(frame_file);
            frame = frame.frame;
            if(dim == 3)
                [lbld, I] = frame.GetPadded(w, ~isempty(refMasked), center(3));
            else
                [lbld, I] = frame.GetPadded(w, ~isempty(refMasked));
            end
            
            if(~isempty(dc))
                I = double(I);
                I = I - mean(I(:)) + dc;
                I(I < 0) = 0;
            end
            dc_out = mean(double(I(:)));
            
            subLabeled = double(lbld(center(2)-w:center(2)+w, center(1)-w:center(1)+w));
            subLabeled(imclearborder(subLabeled ~= 0) == 0) = 0;

            patch = double(I(center(2)-w:center(2)+w, center(1)-w:center(1)+w));

            if(isempty(refMasked))  
                % align according to major axis
                % orientation:
                centerLabel = subLabeled(w+1, w+1);
                mask = (subLabeled == centerLabel);
                if((centerLabel==0) || (max(mask(:))==0)) 
                    %underseg:
                    softMask = mask;
                    status = 0;
                    return;
                end
                
                s = regionprops(mask,'Orientation');
                orientation = cat(1, s.Orientation);
                orientation = orientation(~isnan(orientation));
                deg = 90-orientation(1);
            else
                % nearest non zero label (to patch center):
                centerLabel = -1;
                [~, ind] = bwdist(subLabeled~=0);
                if(max(ind(:)))
                    ind = ind(w+1, w+1);
                    [cY, cX] = ind2sub(size(subLabeled), ind);
                    centerLabel = subLabeled(cY, cX);
                end
                mask = (subLabeled == centerLabel);
                if((centerLabel==0) || (max(mask(:))==0)) 
                    %underseg:
                    softMask = mask;
                    status = 0;
                    return;
                end
                
                if(~isempty(deg))
                    % center cell:
                    s = regionprops(mask,'Centroid', 'Area');
                    subCenter = cat(1, s.Centroid);
                    subArea = cat(1, s.Area);
                    subCenter = subCenter(~isnan(subCenter(:,1)) & ~isnan(subCenter(:,2)) & subArea>1, :);
                    if(isempty(subCenter))
                        softMask = mask;
                        status = 0;
                        return;
                    end
                    if(size(subCenter, 1) > 1)
                        [~, max_ind] = max(subArea);
                        subCenter = subCenter(max_ind, :);
                        warning(['More than one cell with the same label: frame #', num2str(frame.Number), ' label: ', num2str(centerLabel)]);
                    end
                    dt = [w+1, w+1] - subCenter;
                    RA = imref2d(size(patch));
                    patch = imtranslate(patch, RA, dt, 'cubic');
                    mask = imtranslate(mask,RA ,dt, 'cubic');
                    patch(patch < 0) = 0;
                    mask = imclose((mask ~= 0), se);
                end
            end
            
            if(~isempty(deg))
                patch = imrotate(patch, deg, 'bicubic', 'crop');
                mask = imrotate(mask, deg, 'bicubic', 'crop');
                patch(patch < 0) = 0;
                mask = imclose((mask ~= 0), se);
            end
            
            % weights:
            [~,~,~,softMask] = Frame.TransformLabeledImage(mask, 1);
            maxMask = max(softMask(:));
            if(maxMask ~= 0)
                softMask = double(softMask)/max(softMask(:));
            end
            masked = patch .* softMask;
            if(registerFlag && ~isempty(refMasked))
                Rmoving = imref2d(size(masked));
                Rfixed = imref2d(size(refMasked));
                try
                    tform  = imregtform(masked, refMasked, 'rigid', optimizer, metric);

                    patch = imwarp(patch, Rmoving, tform, 'cubic', 'OutputView',Rfixed);
                    mask = imwarp(mask, Rmoving, tform, 'cubic', 'OutputView',Rfixed);

                    patch(patch < 0) = 0;
                    mask = imclose((mask ~= 0), se);
                catch EX
                    disp(['frame #', num2str(frameNum), ' label ', num2str(centerLabel),': ', EX.message]);
                end
                [~,~,~,softMask] = Frame.TransformLabeledImage(mask, 1);
                maxMask = max(softMask(:));
                if(maxMask ~= 0)
                    softMask = softMask/maxMask;
                else
                    status = 0;
                end
            end
        end
        
        function [mother_data, daughter_data] = calcCellPatches(framesPath, curFrameNum,w, center, dim, registerFlag)
            daughter_data = [];
            mother_data = [];
            
            [curPatch, center_weights, center_mask, center_label, status, ref_deg, dc] = Mitodix.calcFrameCellPatches(framesPath, curFrameNum, center, w, dim, [], registerFlag);
            if(~status)
                return;
            end
            
            if(center_label == 0)
                warning(['center label 0 for COM=', num2str(center), ' in t=', num2str(curFrameNum)]);
            end
                
            curMasked = curPatch .* center_weights;
            [prevPatch, ~, prevMask, ~, status] = Mitodix.calcFrameCellPatches(framesPath, curFrameNum-1, center, w, dim, curMasked, 0, ref_deg, dc);
            if(status)
                daughter_data.birth_frame = curFrameNum;
                daughter_data.patches = [prevPatch(:)', curPatch(:)'];
                [~,~,~,wts] = Frame.TransformLabeledImage(center_mask | prevMask, 1);
                [~, ~, ~, daughter_data.features] = wcorr2(curPatch, prevPatch, wts, 1);
                daughter_data.centers = center;
                daughter_data.centers(1:2) = center(1:2) - [w, w];
                daughter_data.daughterLabels = center_label;
                if(center_label == 0)
                    warning(['center label 0 for COM=', num2str(center), ' in t=', num2str(curFrameNum)]);
                end
            end
            
            [nextPatch, ~, nextMask, ~, status] = Mitodix.calcFrameCellPatches(framesPath, curFrameNum+1, center, w, dim, curMasked, registerFlag, ref_deg, dc);
            if(status)
                mother_data.birth_frame = curFrameNum+1;
                mother_data.patches = [curPatch(:)', nextPatch(:)'];
                [~,~,~,center_weights] = Frame.TransformLabeledImage(center_mask | nextMask, 1);
                maxMask = max(center_weights(:));
                if(maxMask ~= 0)
                    center_weights = center_weights/maxMask;
                end
                [~, ~, ~, mother_data.features] = wcorr2(curPatch, nextPatch, center_weights, 1);
                mother_data.centers = center;
                mother_data.centers(1:2) = center(1:2) - [w, w];
                mother_data.motherLabels = center_label;
                mother_data.weights = center_weights(:)';
                %mother_data.motherImages = curPatch;
                %mother_data.motherMasks = center_mask;
            end
        end
        
        % Calculate initial patch frame candidates structure
        function [mother_data, daughter_data] = calcInitialPatchFrameCandidates(framesPath, frameNum, w, registerFlag)
            if(~exist('registerFlag', 'var') || isempty(registerFlag))
                registerFlag = 1;
            end
            curFrame = load(fullfile(framesPath, ['frame', num2str(frameNum), '.mat']));
            curFrame = curFrame.frame;
            centers = curFrame.GetComs(-1);%(w);
            dim = size(centers, 2);
            coms = centers;
            coms(:, 1) = coms(:, 1) + w;
            coms(:, 2) = coms(:, 2) + w;
            cellsNum = size(coms, 1);
            mother_data = [];
            daughter_data = [];
            for c=1:cellsNum
                if(size(coms, 2)==3)
                    com = coms(c,:,:);
                else
                    com = coms(c,:);
                end
                [cell_mother_data, cell_daughter_data] = Mitodix.calcCellPatches(framesPath, frameNum,w, com, dim, registerFlag);
                mother_data = [mother_data; cell_mother_data];
                daughter_data = [daughter_data; cell_daughter_data];
            end
        end
    end
end    

