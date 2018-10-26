classdef Frame < handle
    % This class holds a data about a frame.
    
    properties(Constant)
        Log = SingleInstance.Logger;    % Application logger.
        Configure = Config();           % Application configuration parameters.
    end
    
    properties (Access = private)
        % Cells data:
        %cells;                  % Hashtable of Cell objects.
        centers;                % Cells centers. Nx2 matrix.
        labels;                 % an Nx1 vector of cells' labeles.
        ages;                   % an Nx1 vector of cells' ages [#frames].
        neighbors;              % Mx2 matrix of all neigboring cell pairs (cell labels).
        neighborsDists;         % Mx1 vector of all neigboring cell pairs distances.
        neighborsMidPoints;     % Mx2 matrix of all neigboring cell mid points.
        lengths;                % an Nx1 vector of cells major axis lengths.
        widths;                 % an Nx1 vector of cells semi major axis.
        sizes;                  % Nx1 vector of cells size (in pixels).
        orientations;           % Nx1 vector of cells major axes orientations [deg].
        backGroundDistParams;   % Mean + Var of the background Gaussian noise fit.
    end
    
    properties (SetAccess = protected)
        Number;                 % frame number (T).
        Name;                   % frame file name.
        LabeledName;            % frame segmentation file name.
        Image;                  % Gray scale image of the frame.
        Labeled;                % Labeled segmentation of the frame.
        Contour;                % Labeld contour of the frame.
        BwContour;              % BW contour of the frame.
        Size;                   % Frame size (Y,X).
        SaturationFlag=0;       % Indicate intensity saturation in input gray scale image.
        Z_dim=1;              	% Z dimension (equals 1 for 2D).
        
        % frame statistics:
        CellLength;             % frame cells lengths statistics struct.
        CellWidth;              % frame cells widths statistics struct.
        CellSize;               % frame cells sizes statistics struct.
        CellDrift;              % frame cells drifts statistics struct.
        DistBetweenNeigbors;    % frame cells distances between neigbors statistics struct.
        
        %Image padded with zeros:
%         PaddedImage;            % Gray scale image of the frame padded with zeros.
%         PaddedLabeled;          % Labeled segmentation of the frame padded with zeros.
%         W;                      % Zero padding margin.
    end
    
    properties (Dependent)
        CellsData;              % Cells data structure.
        Center;                 % Frame center (Y,X).
        NumberOfCells;          % Number of cells in frame
        SoftMask;               % Fuzzy segmentation map.
        BW;                     % BW segmentation of the frame.
        Cells;                  % get array of all cells as CellInstance.
        MaxVal;                 % max Image intensity value.
        MinVal;                 % min Image intensity value.
        BgColor;                % Frame mean back ground gray level.
    end
    
    methods
        % Constructor:
        function self = Frame(varargin)
        % Constructor
        %   self = Frame(varargin)
        % Inputs:
        %   First option: input is a struct.
        %   Second option:
        %       <frameNum> = frame number.
        %       <image> = frame's gray scale image.
        %       <labeled> = frame's labeled segmentation.
        %                   If it is BW and not labeled, lebels are 
        %                   calculated according to connected components.
        %       <frameName> =   frame file name.
        %       <COMs> = Hashtable (containers.Map) of cell 
        %                centers [x,y] with cell label as key.
        %       <labeledName> = frame labeled file name.
        %       <rescaleFlag> = 0/1 rescale to uint8 (Default = 1).
        
        % Struct constructor:
            if(isempty(self.Configure.DataName))
                msg = 'Invalid configuration file';
                    self.Log.mlog = {LogType.ERROR ...
                          ,mfilename('class') ...
                          ,[self.Log.Me,msg]};
                    error('Frame:InvalidInput', msg);
            end
            
            if(nargin ==1)
                self.Load(varargin{1});
                return;
            end

        % Regular constructor:
            self.Number = varargin{1};
            self.Name = varargin{4};
            image = varargin{2};
            labeled = varargin{3};
            if((nargin > 5) && ~isempty(varargin{6}))
                self.LabeledName = varargin{6};
                tI = Frame.GetFileNameToken( self.Name, self.Configure.InputImExpression );
                tL = Frame.GetFileNameToken( self.LabeledName, self.Configure.InputSegExpression );
                if(~isequal(tI, tL))
                    msg = ['Frame image file ', self.Name, ' was paired with segmentation file ', self.LabeledName,'.'];
                    self.Log.mlog = {LogType.WARNING ...
                          ,mfilename('class') ...
                          ,[self.Log.Me,msg]};
                    warning('Frame:InvalidInput', msg);
                end
            end
            
            if((nargin > 6) && ~isempty(varargin{7}))
                rescaleFlag = varargin{7};
            else
                rescaleFlag = 1;
            end
            
            if((nargin > 7) && ~isempty(varargin{8}))
                minCellLengthToConsider = varargin{8};
            else
                minCellLengthToConsider = self.Configure.MinCellLength;
            end
            
			% Do not consider small regions as cells:
            if(isempty(minCellLengthToConsider))
                minCellLengthToConsider = 4;
            end
			
            %self.cells = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
			if(length(size(image))==3)
				self.Z_dim = size(image, 3);
			else
				self.Z_dim = 1;
			end
			Y = size(image, 1);
			X = size(image, 2);
            self.Size = [Y, X];
			if(self.Z_dim > 1)
				labeled = labeled(1:Y, 1:X, 1:self.Z_dim);
			else
				labeled = labeled(1:Y, 1:X);
            end
            lblVec = reshape(labeled,[],1);
            % Remove empty z slices:
            if(self.Z_dim > 1)
                for zz=1:self.Z_dim
                    lbld = labeled(:, :, zz);
                    if any(lbld(:))
                        break;
                    end
                end
                min_z = max(1, zz-1);
                for zz=self.Z_dim:-1:1
                    lbld = labeled(:, :, zz);
                    if any(lbld(:))
                        break;
                    end
                end
                max_z = min(self.Z_dim, zz+1);
                image = image(:,:,min_z:max_z);
                labeled = labeled(:,:,min_z:max_z);
                self.Z_dim = size(image, 3);
                lblVec = reshape(labeled,[],1);
            end
            
            % input segmentation is binary or labeled:
            isLabeled = ((sum(lblVec ~= 0) - sum(lblVec == max(lblVec))) > 0);
            
            % Watershed:
            if(~isLabeled)
                [~, split] = SplitBigClusters(labeled);
            else
                split = zeros(size(labeled));
                %[~, split] = WaterShed(labeled);
            end
            labeled(split~=0) = 0;
            relabelFlag = any(split(:));

            % connected componenets:
            if((self.Z_dim > 1) || (~isLabeled))
                CC = bwconncomp(labeled~=0);
                lbld = labelmatrix(CC);
            else
                CC = labeled;
                lbld = labeled;
            end
            
            if(self.Z_dim > 1)
                s = regionprops(CC,'Area');
                maLengths = ((3/(4*pi))*cat(1, s.Area)).^(1/3); % radi according to sphere volume
            else
                s = regionprops(CC,'MajorAxisLength');
                maLengths = cat(1, s.MajorAxisLength);
            end
            
            % Check if major axis min length filter is too big:
            maxLengths = max(maLengths);
            if(minCellLengthToConsider > maxLengths)
                self.Log.mlog = {LogType.WARNING ...
                          ,mfilename('class') ...
                          ,[self.Log.Me,' min cell diameter filter is too big in frame ', num2str(self.Number),' (minCellLengthToConsider = ',num2str(minCellLengthToConsider),'[pixels] while maximal major axis in this frame is ',num2str(maxLengths),').']};
            end
            
            % remove noise segmentation:
            idx = find(maLengths <= max(1,minCellLengthToConsider));%min(4,minCellLengthToConsider));
            if(~isempty(idx))
                for ii=1:length(idx)
                    labeled(lbld == idx(ii)) = 0;
                    lbld(lbld == idx(ii)) = 0;
                end
            end
            
            if(relabelFlag || ~isLabeled)
                labeled = lbld;
            end
            
            self.Labeled = labeled;
			lblVec = reshape(self.Labeled, [], 1);
			
            % Mask border cells from cells list:
			innerLabels = labeled;			
% 			innerLabels(imclearborder(labeled ~= 0) == 0) = 0;
			
			% Get connected components properties:
            s = regionprops(innerLabels, 'Centroid','Area');
            self.sizes = cat(1, s.Area);
            
            if(self.Z_dim == 1)
                s2 = regionprops(innerLabels, 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
                self.lengths = cat(1, s2.MajorAxisLength);
                self.widths = cat(1, s2.MinorAxisLength);
                self.orientations = cat(1, s2.Orientation);
            else
                self.lengths = ((3/(4*pi))*cat(1, s.Area)).^(1/3);
                self.widths = self.lengths;
                self.orientations = zeros(size(self.lengths));
            end

            idx = find((self.sizes>0) & (self.lengths > minCellLengthToConsider) & (self.sizes > 2*minCellLengthToConsider));
            self.centers = round(cat(1, s(idx).Centroid));
            self.widths = self.widths(idx);
            self.orientations = self.orientations(idx);
            self.sizes = self.sizes(idx);
            self.lengths = self.lengths(idx);
            self.labels = idx;
            
            bgMean = mean(image(lblVec==0));
            N = length(idx);
            max_label = max(self.Labeled(:));
            labeled = self.Labeled;
            self.Labeled = 0*labeled;
            for n=1:N
                label = self.labels(n);

                % over seg:
                if(mean(image(lblVec == label)) < bgMean)
                    %labeled(labeled == label) = 0;
                    self.labels(n) = -1;
                    continue;
                end
                
				%if(self.Z_dim == 1)
                % make sure only one blob has this label:
                sub_l = labelmatrix(bwconncomp(labeled == label));
                if(max(sub_l(:)) > 1)
                    if(self.Z_dim == 1)
                        sub_s = regionprops(sub_l, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
                        sub_areas = cat(1, sub_s.Area);
                        sub_lengths = cat(1, sub_s.MajorAxisLength);
                        sub_widths = cat(1, sub_s.MinorAxisLength);
                        sub_orientations = cat(1, sub_s.Orientation);
                    else
                        sub_s = regionprops(sub_l, 'Area', 'Centroid');
                        sub_areas = cat(1, sub_s.Area);
                        sub_lengths = ((3/(4*pi))*cat(1, sub_s.Area)).^(1/3);
                        sub_widths = sub_lengths;
                        sub_orientations = zeros(size(sub_lengths));
                    end
                    
                    sub_centers = round(cat(1, sub_s.Centroid));
                    sub_idx = find((sub_lengths > max(1,minCellLengthToConsider)) & (sub_areas > 2*max(1,minCellLengthToConsider)));
                    small_idx = find((sub_lengths <= max(1,minCellLengthToConsider)) | (sub_areas <= 2*max(1,minCellLengthToConsider)));
                    if(isempty(sub_idx))
                        self.labels(n) = -1;
                        %self.Labeled(sub_l ~= 0) = 0;
                    else
%                         % Remove small parts:
%                         for area_k=1:length(small_idx)
%                             self.Labeled(sub_l ~= 0 & sub_l == small_idx(area_k)) = 0;
%                         end
                        
                        % Set largest part as a cell:
                        sub_areas = sub_areas(sub_idx);
                        sub_lengths = sub_lengths(sub_idx);
                        sub_widths = sub_widths(sub_idx);
                        sub_orientations = sub_orientations(sub_idx);
                        sub_centers = sub_centers(sub_idx, :);
                        
                        % Set each part as a different cell:
                        [~, areas_idx] = sort(sub_areas, 'descend');
                        self.widths(n) = sub_widths(areas_idx(1));
                        self.orientations(n) = sub_orientations(areas_idx(1));
                        self.sizes(n) = sub_areas(areas_idx(1));
                        self.lengths(n) = sub_lengths(areas_idx(1));
                        self.centers(n,:) = sub_centers(areas_idx(1), :);
                        self.Labeled(sub_l == areas_idx(1)) = label;
                        self.labels(n) = label;
                        
                        for area_k=2:length(areas_idx)
                            max_label = max_label + 1;
                            self.Labeled(sub_l == areas_idx(area_k)) = max_label;
                            self.labels(end+1) = max_label;
                            self.widths(end+1) = sub_widths(areas_idx(area_k));
                            self.orientations(end+1) = sub_orientations(areas_idx(area_k));
                            self.sizes(end+1) = sub_areas(areas_idx(area_k));
                            self.lengths(end+1) = sub_lengths(areas_idx(area_k));
                            self.centers(end+1,:) = sub_centers(areas_idx(area_k), :);
                        end
                    end
                else
                    self.labels(n) = label;
                    self.Labeled(labeled == label) = label;
                end
            end
            
            
            % remove over seg:
            ind = (self.labels ~= -1);
            self.widths = self.widths(ind);
            self.sizes = self.sizes(ind);
            self.orientations = self.orientations(ind);
            self.lengths = self.lengths(ind);
            self.labels = self.labels(ind);
            self.centers = self.centers(ind, :);
            
            %% Use better centers if avilable:
            %if(isLabeled && (nargin > 4) && ~isempty(varargin{5}))
            %    COMs = varargin{5};
            %    lbls = keys(COMs);
            %    for k=1:length(lbls)
            %        key = lbls(k);                    
            %        lbl = str2num(key);
			%
            %        % if label exists in labels vector, replace the
            %        % center:
            %        idx = (self.labels == lbl);
            %        if(count(idx))
            %            self.centers(idx) = COMs(key);
            %        end
            %    end
            %end

            % Calculate statistics:
            self.CellLength = Statistics(self.lengths);
            self.CellWidth = Statistics(self.widths);
            self.CellSize = Statistics(self.sizes);
            
            % check saturation:
            [h, v] = hist(double(image(image(:)>0)), 256);
            self.SaturationFlag = (h(end) > mean(h(end-5:end-1))*100);
            
            self.Image = image;
            
            % Updated background estimation:
%             bgPixels = image(self.Labeled == 0);
%             bgMu = mean(bgPixels);
%             bgSigma = var(bgPixels);
            self.backGroundDistParams = []; 

            % Find negibors:
            self.setNeighbors();
            self.DistBetweenNeigbors = Statistics(self.neighborsDists);
            
            % padd the images with zero margins:
%             self.W = [];%round(self.CellLength.Max + self.CellLength.STD);
%             self.PaddedImage = [];
%             self.PaddedLabeled = [];

%             if(rescaleFlag)
%                 self.Image = uint8(255*self.Image);
%                 self.PaddedImage = uint8(255*self.PaddedImage);
% %                 self.Labeled = uint16(self.Labeled);
% %                 self.PaddedLabeled = uint16(self.PaddedLabeled);
%             end
            self.ages = self.Configure.CellCycle*ones(length(self.labels),1);
        end

        % Get CellsData struct
        function CellsData = get.CellsData(self)
            CellsData = [];
            CellsData.Centers = self.centers;
            CellsData.Labels = self.labels;               
            CellsData.Ages = self.ages;                 
            CellsData.Neighbors = self.neighbors;              
            CellsData.NeighborsDists = self.neighborsDists;
            CellsData.NeighborsMidPoints = self.neighborsMidPoints; 
            CellsData.Lengths = self.lengths;               
            CellsData.Widths = self.widths;                
            CellsData.Sizes = self.sizes;
            CellsData.Orientations = self.orientations; 
        end

        % Get maximal Image intensity value
        function MaxVal = get.MaxVal(self)
            MaxVal = max(self.Image(:));
        end
        
        % Get minimal Image intensity value
        function MinVal = get.MinVal(self)
            MinVal = min(self.Image(:));
        end
        
        % Get mean background value
        function BgColor = get.BgColor(self)
            BgColor = mean(self.Image(self.Labeled(:)>0));
        end
        
        % Get frame center
        function Center = get.Center(self)
            Center = Frame.GetFrameCenter(self.Size(2), self.Size(1));
        end
        
        % Get number of cells in frame
        function NumberOfCells = get.NumberOfCells(self)
            NumberOfCells = length(self.labels);
        end
        
        % Get BW
        function BW = get.BW(self)
            BW = (self.Labeled > 0);
        end
        
        % Get Fuzzy segmentation map
        function SoftMask = get.SoftMask(self)
            [~, ~, ~, SoftMask] = Frame.TransformLabeledImage(self.Labeled);
        end
        
        % Get all cells
        function Cells = get.Cells(self)
            Cells = [];
            for l=1:length(self.labels)
                Cells = [Cells, self.GetCell(self.labels(l))];
            end
        end
        
        % Get padded image and labeled image
        function [PaddedLabeled, PaddedImage] = GetPadded(self, w, removeBorders, z_idx)
            if(~exist('removeBorders', 'var') || isempty(removeBorders))
                removeBorders = 1;
            end
            
            if(~exist('z_idx', 'var') || isempty(z_idx) || (self.Z_dim == 1))
                im = self.Image;
                lbld = self.Labeled;
            else
                im = self.Image(:, :, z_idx);
                lbld = self.Labeled(:, :, z_idx);
            end
            bgColor = mean(im(lbld(:) == 0));
            
            if(removeBorders)
                innerLabels = imclearborder(lbld ~= 0);               
                im((lbld ~= 0) & (innerLabels == 0)) = bgColor;
                lbld(innerLabels == 0) = 0;
            end

            PaddedImage = padarray(im, [w, w], bgColor, 'both');
            PaddedLabeled = padarray(lbld, [w, w], 0, 'both');
        end
        
        % get inner cells centers
        coms = GetComs(self, margin)
        
        % load struct data
        Load(self, sobj)
        
        % Perform subclass save operations
        sobj = saveobj(self)
        
        % Set cells drifts
        updateFrame = SetCellDrift(self, prevFrameCenters)
        
        % Update frame data
        Update(self, varargin)
        
        % Get cell object
        cell = GetCell(self, cellLabel)
        
        % Get cell data
        [center, length, width, size, orientation] = GetCellData(self, cellLabel)
        
        % Get ages of input labeles cells
        ages = GetAges(self, lbls)
        
        % Set cells ages
        SetAges(self, ages)

        % set frame image back ground color (gray level of bg pixels)
        SetBgColor(self, bgColor)
        
        % Get data of neigbors cells in a search radius
        [neighborsLabels, neighborsMidPoints, neighborsDists] = GetNeigbors(self, R, cellLabels)
        
        % Get cell labels in a search radius
        [labels, dists] = GetCellLabelsInR(self, ref, R, k)
        
        % Check if input label is a border cell
        isOnBorder = IsLabelOnBorder(self,cellLabel)
        
        % Estimate cell shape on the border of the frame
        [isOnBorder, img, lbld] = EstimateBorderCell(self, cellLabel)
        
        % Get closes cell label
        [labels, idx] = GetClosestCellLabel(self, newCOMs, newLabels)
        
        % Get cells with labeles bigger than the maximal in prev frame.
        [newCenters, newLabels] = GetNewerCellsCenters(self, maxPreviousLabelValue)
        
        % Get possible mothers
        [motherLabeles, foundFlags, R] = GetPossibleMothers(self, daughtersMidPoints, R, anchorMothers)
        
        % Get search radius according to cell drift statistics
        R = GetDriftRadius( self )
        
        % Get search radius according to cell distance between neigbhors statistics
        R = GetDistRadius( self )
        
        % Rescale Gray level:
        RescaleGrayLevel( self, factor, uintFlag, offset )
        
        % Get anchor cells for mitodix search:
        [ labels, scores ] = GetAnchorCells( self, prevFrame, R, Rm )
        
        % Get closest pixel label
        function label = ClosestLabel(self, com)
            label = Frame.GetClosestLabel(self.Labeled, com);
        end
        
        % Get patch around midPoint containing only cells of requested 2 labels where the cells are one above the other.
        [ im, weights, mask, cellsImages, coms, d_deg ] = GetCellsPatch( self, labels, pSize, alignFlag, divideFlag )
    end
    
    methods (Access = private)
        % Re-calculate cell data
        setCellData(self, cellLabel)
        
        % Set neigbors graph in frame
        setNeighbors( self )
        
        % Add neigbor to neigbors list only if was not added
        addNeighbor(self, id1, id2, c1, c2)
    end
    
    methods (Static)
        % Load object from struct
        function obj = loadobj(objStruct)
            obj = Frame(objStruct);
        end
        
        % Get frame center
        frameCenterPoint = GetFrameCenter(X, Y)
        
        % Get closest pixel label
        label = GetClosestLabel(labeled, com)
        
        % Find element in double sorted array
        [index, found] = FindDoubleSortedElementPosition(array, val)
        
        % Calculate image transformation
        [imContour, imSdf, imPhi, imSoftMask, biggerMask] = TransformLabeledImage(labeled, epsilon, mask, softThr, mskThr)
        
        % Get image bounding box
        bbox = GetBBox(image, COMyx)
        
        % Re-calculate image segmentation
        seg = CalcSeg(image, labeled)
        
        % Save frames patch images
        SaveImages(frames, outputPath, com, R)
        
        % Re-calculate segmentation to separate connected cells
        bw = SeparateCellsSeg(seg)
        
        % Perform segmentation re-calculation on an images batch
        SeparateCellsSegBatch(mtdx, startSegFrame, endSegFrame)
        
        % Calculate mask (round edges)
        mask = CalcMask(BW)
        
        % Get file names token according to input expression
        [t, fileNames] = GetFileNameToken( fileNames, expr )
        
        % Calculate weighted similarity between 2 patches.
        sim = PatchSim( im1, mask1, im2, mask2, registerFlag )
        
        % Angle between the line that connects the two COMs and the frame's Y axis
        % under the assumption that daughter1 has smaller Y coordinates than daughter2.
        orientation = CalcOrientation(dA, dB)
        
        % Read intensity image and lableled image from file paths (2D or 3D).
        [ image, labeled, Z_dim ] = Read2D3D( im_path, seg_path, voxel_idx, rescaleFlag, downscaleFactor )
        
        % Create Disk segmentations.
        DisksSeg( dataName, outputPath, radi )
    end
end

