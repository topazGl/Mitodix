function tracks = Tracking(inDir, outDir, full_file_names, FOI_E, useFillHoles, minSegmAreaPx, Rm, minOvlp, mitosisCandidates)
%%%
%	This function generates the tracking presented in:
%	Gilad et al, Fully Unsupervised Symmetry-Based Mitosis Detection in Time-Lapse Cell Microscopy, 2018.
%   and is based on the followin code of the FR-Ro-GE group of the cell tracking challange:
%   http://www.celltrackingchallenge.net/participants/FR-Ro-GE/
%%%
%   Example:
% inDir = '..\Data\Outputs\tracking\Fluo-N2DL-HeLa\01\in_seg\';   % input segmentations.
% outDir = '..\Data\Outputs\tracking\Fluo-N2DL-HeLa\01\tracking'; % output tracking.
% fileName = dir(fullfile(inDir, '*.tif'));
% [tokens, full_file_names] = Frame.GetFileNameToken({fileName.name}, 't([0-9]+)_Simple Segmentation.tif');
% expr = 't([0-9]+)*.tif';
% Rm = 50; minOvlp = 0.7;
% tracks = Tracking(inDir, outDir, full_file_names, [], [], [], Rm, minOvlp);
%% output video:
% inDirIntensity = '..\Data\Inputs\ISBI\2014\Fluo-N2DL-HeLa\01\'; % input gray-level images (for video).
% TrackingVideo(inDirIntensity, outDir, expr, full_file_names, outDir, tracks, 30, 0, 0);
% TrackingVideo(inDirIntensity, outDir, expr, full_file_names, fullfile(outDir, 'divisions'), tracks, 30, 0, 1);
%
%% To evalutae mitosis detection performance:
% %create res_track_4c.txt!
% performance = evalTRA(outDir, '..\Data\Inputs\ISBI\2014\Fluo-N2DL-HeLa\01_GT\GT\');
% save(fullfile(outDir, 'performance.mat'), 'performance');

    if(~exist('FOI_E', 'var') || isempty(FOI_E))
        FOI_E = 0;
    end
    
    if(~exist('useFillHoles', 'var'))
        useFillHoles = [];
    end
    
    if(~exist('minSegmAreaPx', 'var'))
        minSegmAreaPx = [];
    end
    
    if(~exist('Rm', 'var'))
        Rm = [];
    end
    
    if(~exist('minOvlp', 'var'))
        minOvlp = [];
    end
    
    if(~exist('mitosisCandidates', 'var'))
        mitosisCandidates = [];
    end
    
    if(~exist(outDir, 'dir'))
        mkdir(outDir);
    end
    
    % Initialize tracks
    tracks = [];
    maskFOI = [];
    L2_prev = [];
    bg_lbl = [];
    binary_input_flag = 0;
    nFrames = length(full_file_names);
	for frameIdx=1:nFrames
        fn = full_file_names{frameIdx};
        L2 = imread(fullfile(inDir, fn));
        [L2, bg_lbl] = preprocessBinarySegmentation(L2, useFillHoles, minSegmAreaPx, bg_lbl); 
        [L2, tracks, maskFOI, L2_FOI] = propogateTracking(frameIdx, L2, L2_prev, tracks, FOI_E, maskFOI, Rm, minOvlp, mitosisCandidates);
        L2_prev = L2;
        
      %
      % Write segmentation mask (propagated labels, uint16)
      %
      fnOut = fullfile(outDir, fn);
      imwrite(uint16(L2_FOI), fnOut,'WriteMode','overwrite');
	end

    %
    % Store tracking result
    %
    fnOut = fullfile(outDir, 'res_track.txt');
    fprintf('\twrite %s\n', fnOut);

    fileID = fopen(fnOut, 'w');
    try
        nTracks = length(tracks);
        for i = 1:nTracks
          if (tracks(i, 1) ~= -1)
              fprintf(fileID, '%d %d %d %d %.4f %.4f\n', tracks(i,1), tracks(i,2), tracks(i,3), tracks(i,4), tracks(i,5), tracks(i,6));
    %           fprintf('%d %d %d %d %.4f\n', tracks(i,1), tracks(i,2), tracks(i,3), tracks(i,4), tracks(i,5), tracks(i,6));
          else
              fprintf('track %d is empty\n', i);
          end
        end
        fclose(fileID);
    catch EX
        fclose(fileID);
        rethrow EX;
    end
end

function [L2, org_bg_lbl] = preprocessBinarySegmentation(BW, useFillHoles, minSegmAreaPx, org_bg_lbl)
  if(~exist('useFillHoles', 'var') || isempty(useFillHoles))
      useFillHoles = 0;
  end
  if(~exist('minSegmAreaPx', 'var') || isempty(minSegmAreaPx))
    minSegmAreaPx = 100;
  end
  
  %
  % make sure BG is '0':
  %
  if(~exist('bg_lbl', 'var') || isempty(org_bg_lbl))
      org_bg_lbl = 0;
  end
  BW = (BW ~= org_bg_lbl);
  
  %
  % Fill holes
  %
  if (useFillHoles==1)
    BW = imfill(BW,'holes');
  end
  CC = bwconncomp(BW);
  
  %
  % Remove small segments
  %
  if (minSegmAreaPx > 0)
    numPixels = cellfun(@numel,CC.PixelIdxList);
    CC2 = CC;
    CC2.PixelIdxList = CC2.PixelIdxList(numPixels>minSegmAreaPx);
    CC2.NumObjects = length(CC2.PixelIdxList);   
    CC = CC2;
  end  
  
  L2 = uint16(labelmatrix(CC));
end

function [L2, tracks, maskFOI, L2_FOI] = propogateTracking(frameIdx, L2, L2_prev, tracks, FOI_E, maskFOI, Rm, minOvlp, mitosisCandidates)
    if(~exist('Rm', 'var') || isempty(Rm))
        Rm = 0.1*min([size(L2, 1), size(L2, 2)]);
    end
    
    if(~exist('minOvlp', 'var') || isempty(minOvlp))
        minOvlp = 0;
    end
    
    if(~exist('mitosisCandidates', 'var'))
        mitosisCandidates = [];
    end
  %
  % Remove objects that are completely outside the FOI
  %
  L2_FOI = L2;
  if (FOI_E > 0)
      if (frameIdx==1)
          % Initialize FOI mask
          maskFOI = false(size(L2));
          maskFOI(FOI_E+1:end-FOI_E, FOI_E+1:end-FOI_E) = true;
      end
      rp2 = regionprops(L2);
      L2_tmp = L2;
      L2_tmp(~maskFOI) = 0;
      rp3 = regionprops( L2_tmp);

      for i = 1:length(rp2)
        if (i<=length(rp3))
            if (rp3(i).Area==0)
              L2_FOI(L2_FOI==i) = 0;
            end
        else
            L2_FOI(L2_FOI==i) = 0;
        end
      end
  end
  
  %
  % Tracking / Label propagation
  %
  curr_labels = sort(unique(L2(:)));
  curr_labels = curr_labels(curr_labels>0);
  curr_nObjects = length(curr_labels);
  if (frameIdx==1)
    % Initialize objects and labels  
    %assert(curr_nObjects>0,'assert failed: curr_nObjects>0; No objects at initialization step.');
    
    % Initialize tracks
    nTracks = curr_nObjects;
    tracks = zeros(nTracks, 6, 'double');
    for i = 1:nTracks
        if (sum(sum(L2_FOI(L2==i)>0)) > 0)
            tracks(i, :) = [i, frameIdx-1,frameIdx-1, 0, double(1), frameIdx-1];
        else
            tracks(i, 1) = double(-1);
        end
    end
  else
    % Propagate object labels from previous segmentation mask
    lblMapping = zeros(curr_nObjects, 1);
    
    prev_labels = sort(unique(L2_prev(:)));
    prev_labels = prev_labels(prev_labels>0);
    nObjects = length(prev_labels);
  
    if (nObjects>=0 && curr_nObjects>0)
    
        if (nObjects>0)
            %
            % compute intersection over union
            %
            currLin = double(L2(:))+1;
            prevLin = double(L2_prev(:))+1;
            intersec = sparse(prevLin, currLin, ones(size(prevLin)), ...
                                        max(prevLin), max(currLin));
            intersec = full(intersec);    
            areaPrev = sum(intersec,2);
            areaCurr = sum(intersec,1);
            
            intersec = intersec(2:end,2:end);
            areaPrev = areaPrev(2:end);
            areaCurr = areaCurr(2:end);
            
            lblMatrix = zeros(nObjects,curr_nObjects);
            
            idx = find(intersec>0);
            [idx_i,idx_j] = ind2sub(size(intersec), idx);
            for k = 1:length(idx_i)
                i = idx_i(k);
                j = idx_j(k);
                areaUnion = areaPrev(i) + areaCurr(j) - intersec(i,j);
                lblMatrix(find(prev_labels==i, 1, 'first'),find(curr_labels==j, 1, 'first')) = intersec(i,j)/areaUnion;             
            end
            
            % Enforces one-to-one mapping

            % Maximum intersection over union
            % src-objs. > curr-objs.
            [~,max_idx] = max(lblMatrix,[],2);
            linidx = sub2ind(size(lblMatrix),1:size(lblMatrix,1),max_idx');
            tmp_lblMatrix = zeros(size(lblMatrix));
            tmp_lblMatrix(linidx) = lblMatrix(linidx);
            lblMatrix = tmp_lblMatrix;

            % Maximum intersection over union
            % curr-objs.
            [max_val,max_idx] = max(lblMatrix,[],1);
            lblMapping = prev_labels(max_idx);

        else % no cells in previous frame
            max_val = zeros(nObjects, 1);
        end
        
        lbl_likelihood = max_val;
        
        % no overlap between current cell to prev frame:
        fullMapping = lblMapping;
        fullMapping(max_val==0) = 0;
        
        nAddObj = sum(max_val<=minOvlp);
        if(nAddObj)           
            lblMapping(max_val<=minOvlp) = 1:nAddObj;
            lblMapping(max_val<=minOvlp) = lblMapping(max_val<=minOvlp) + max(tracks(:,1));
        end
        
        s = regionprops(L2_prev, 'Centroid');
        prev_centers = cat(1, s.Centroid);
        
        % Propagate matched labels
        tmp_L2 = zeros(size(L2));
        parents = []; no_parents = []; cands = [];
        for jj=1:length(curr_labels)
            % j: curr label
            % lblMapping(jj): propagated label
            tmp_lbl = curr_labels(jj);
            cell_lbl = lblMapping(jj);
            alternative_lbl = fullMapping(jj);
            tmp_L2(L2==tmp_lbl) = cell_lbl;
            %inner_cell = sum(inner_L2(:)==j)>0;

            track_ind = find(tracks(:,1)==cell_lbl, 1, 'first');
            if(isempty(track_ind))
                track_ind = 0;
            end
            
            if (sum(sum(L2_FOI(L2==tmp_lbl)>0)) > 0)
                % Update tracks
                if(track_ind)
                    % Update track:
                	track_entry = tracks(track_ind,:);
                    track_entry(3) = frameIdx-1;
%                     min_ovlp_frame = [frameIdx-1, track_entry(6)];
                    track_entry(5) = lbl_likelihood(jj);
                    track_entry(6) = frameIdx-1;
%                     [track_entry(5), min_ovlp_idx] = min([lbl_likelihood(j), track_entry(5)]);
%                     track_entry(6) = min_ovlp_frame(min_ovlp_idx);
                    tracks(track_ind,:) = track_entry;
                else%if(inner_cell)
                    % initiate a new track:
                    parent = 0;
                    cand_idx = 0;
                    s = regionprops(tmp_L2==cell_lbl, 'Centroid', 'Area');
                    center = cat(1, s.Centroid);
                    
                    if(~isempty(mitosisCandidates))
                        % set mitosis candidate cell as parent:
                        frame_ind = find(cat(1, mitosisCandidates.FrameNum) == frameIdx);
                        frame_cands = mitosisCandidates(frame_ind);
                        if(~isempty(frame_cands))
                            d1_coms = cat(1, frame_cands.Daughter1COM);
                            d2_coms = cat(1, frame_cands.Daughter2COM);
                            dist1 = (center(1) - d1_coms(:, 1)).^2 + (center(2) - d1_coms(:, 2)).^2;
                            dist2 = (center(1) - d2_coms(:, 1)).^2 + (center(2) - d2_coms(:, 2)).^2;
                            dists = min([dist1, dist2],[],2);
                            [D, cand_idx] = min(dists);
                            area = cat(1, s.Area);
                            if(D <= 0.5*sqrt(double(area(1))))
                                parent = L2_prev(frame_cands(cand_idx).COM(2), frame_cands(cand_idx).COM(1));
                                cand_idx = frame_ind(cand_idx);
                            end
                        end
                    end
                    
                    if(cand_idx == 0)
                        % set nearest cell as parent:
                        if((alternative_lbl ~= cell_lbl) && (alternative_lbl ~= 0))
                            parent = alternative_lbl;
                        elseif(isempty(mitosisCandidates))
                            dist = (prev_centers(prev_labels,1)-center(1)).^2 + (prev_centers(prev_labels,2)-center(2)).^2;
                            [D, parent_ind] = min(dist);
                            if(D <= Rm^2)
                                parent = prev_labels(parent_ind);
                            end
                        end
                    end

                    if(parent ~= 0)
                        parents = [parents; parent];
                        cands = [cands; cand_idx];
                    else
                        no_parents = [no_parents, cell_lbl];  
                    end
                        
                    track_ind = size(tracks, 1) + 1;
                    track_entry = [cell_lbl, frameIdx-1, frameIdx-1, parent, double(1), frameIdx-1];
                    tracks(track_ind,:) = track_entry;
                end
            end
        end
        parents = unique(parents);
        new_flag = zeros(length(parents), 1);
        for jj=1:length(parents)
            parent = parents(jj);
            track_ind = find(tracks(:,1)==parent, 1, 'first');
            
            daughters_count = sum(tracks(:,4) == parent);
            if(~daughters_count)
                continue;
            end
            
            if(tracks(track_ind,3) == frameIdx-1)
                if(isempty(mitosisCandidates))
                    tracks(tracks(:,4) == parent, 4) = 0;
                else
                    % terminate parent track:
                    tracks(track_ind,3) = frameIdx-2;

                    % initiate a new one:
                    new_lbl =  max(tracks(:,1))+1;
                    tracks(size(tracks, 1) + 1, :) = [new_lbl, frameIdx-1, frameIdx-1, parent, double(1), frameIdx-1];
                    tmp_L2(tmp_L2==parent) = new_lbl;
                    new_flag(jj) = 1;
                end
            end
        end
        
        for jj=1:length(parents)
            parent = parents(jj);
            track_ind = find(tracks(:,1)==parent, 1, 'first');
            daughter_tracks = find(tracks(:,4) == parent);
            daughters_count = length(daughter_tracks);
            if(~daughters_count)
                continue;
            end
            daughters_count = daughters_count - new_flag(jj);
%             else
                % check if has two daughters:
            daughter_tracks = find(tracks(:,4) == parent);
            if(cands(jj) ~= 0)
                if(daughters_count~=2)
                %Use candidate cells
                    cand = mitosisCandidates(cands(jj));
                    lbl1 = tmp_L2(cand.Daughter1COM(2), cand.Daughter1COM(1));
                    if(lbl1~=0)
                        track_d1 = find(tracks(:,1)==lbl1, 1, 'first');
                        tracks(track_d1,:) = [lbl1, frameIdx-1, frameIdx-1, parent, double(1), frameIdx-1];
                    end
                    lbl2 = tmp_L2(cand.Daughter2COM(2), cand.Daughter2COM(1));
                    if(lbl2~=0)
                        track_d2 = find(tracks(:,1)==lbl2, 1, 'first');
                        tracks(track_d2,:) = [lbl2, frameIdx-1, frameIdx-1, parent, double(1), frameIdx-1];
                    end
                    tracks(daughter_tracks(daughter_tracks~=track_d1 & daughter_tracks~=track_d2), 4) = 0;
                end
            end
            
            % Remove empty tracks:
            daughter_tracks = find(tracks(:,4) == parent);
            daughters_count = length(daughter_tracks);
            for dd=1:length(daughter_tracks)
                if(~any(tmp_L2(:)==tracks(daughter_tracks(dd), 1)))
                    tracks(daughter_tracks(dd), :) = [-1 0 0 0 0 0];
                    daughters_count = daughters_count-1;
                end
            end
            
            daughter_tracks = find(tracks(:,4) == parent);
            daughters_count = length(daughter_tracks);
            if(daughters_count==1)
                % continue original track:
                daughter_lbl = tracks(daughter_tracks, 1);
                tmp_L2(tmp_L2==daughter_lbl) = parent;
                tracks(track_ind, 3) = frameIdx-1;
                tracks(track_ind, 5) = 0;
                tracks(track_ind, 6) = frameIdx-1;
                tracks(daughter_tracks, :) = [-1, 0, 0, 0, 0, 0];
            elseif(daughters_count>2)
                %prefer inner cells as daughters:
                tmp_L2_inner = imclearborder(tmp_L2);
                for dd=1:length(daughter_tracks)
                    if(~any(tmp_L2_inner(:)==tracks(daughter_tracks(dd), 1)))
                        tracks(daughter_tracks(dd), 4) = 0;
                        daughters_count = daughters_count-1;
                    end
                    if(daughters_count == 2)
                        break;
                    end
                end
                if(daughters_count > 2)
                    daughter_tracks = find(tracks(:,4) == parent);
                    daughters_count = length(daughter_tracks);
                    daughter_centers = zeros(daughters_count, 2);
                    for dd=1:daughters_count
                        daughter_lbl = tracks(daughter_tracks(dd), 1);
                        s = regionprops(tmp_L2==daughter_lbl, 'Centroid');
                        daughter_centers(dd,:) = cat(1, s.Centroid);
                    end
                    % set nearest pair as daughters:
                    [Idx,D] = knnsearch(daughter_centers, daughter_centers, 'K', 2);
                    Idx = Idx(:,end); D = D(:, end);
                    [~, d1] = min(D);
                    d2 = Idx(d1);
                    
                    for dd=1:daughters_count
                        if((dd ~= d1) && (dd ~= d2))
                            tracks(daughter_tracks(dd), 4) = 0;
                        end
                    end
                end
            end
            
            daughter_tracks = find(tracks(:,4) == parent);
            daughters_count = length(daughter_tracks);
            if(daughters_count>0 && daughters_count ~= 2)
                warning([num2str(daughters_count), ' cells accosiated to parent cell ', num2str(parent)])
            end
%             end
        end
        
        new_tracks = find(tracks(:,2) == frameIdx-1);
        for dd=1:length(new_tracks)
            if(~any(tmp_L2(:)==tracks(new_tracks(dd), 1)))
                tracks(new_tracks(dd), :) = [-1, 0, 0, 0, 0, 0];
            end
        end
        
        % pair new and terminated non-mitotic tracks:
        new_tracks = find(tracks(:,2) == frameIdx-1 & tracks(:,4) == 0);
        
        if(~isempty(new_tracks))
            terminated_tracks = find(tracks(:,3) == frameIdx-2);
            terminated_count = length(terminated_tracks);
            paired = zeros(terminated_count, 1);
            for dd=1:terminated_count
                if(any(tracks(:, 4)==tracks(terminated_tracks(dd), 1)))
                    paired(dd) = 1;
                end
            end
            terminated_tracks = terminated_tracks(~paired);
            if(~isempty(new_tracks))
                terminated_centers = zeros(length(terminated_tracks), 2);
                for dd=1:length(terminated_tracks)
                    s = regionprops(L2_prev==tracks(terminated_tracks(dd), 1), 'Centroid');
                    center = cat(1, s.Centroid);
                    if(~isempty(center))
                        terminated_centers(dd,:) = center;
                    end
                end
                new_centers = zeros(length(new_tracks), 2);
                for dd=1:length(new_tracks)
                    s = regionprops(tmp_L2==tracks(new_tracks(dd), 1), 'Centroid');
                    center = cat(1, s.Centroid);
                    if(~isempty(center))
                        new_centers(dd,:) = center;
                    end
                end
                [ind_new,D] = knnsearch(new_centers, terminated_centers);
                [~, idx_prev] = sort(D);
                terminated_tracks = terminated_tracks(idx_prev);
                ind_new = ind_new(idx_prev);
                new_tracks = new_tracks(ind_new);
                for dd=1:length(new_tracks)
                    new_track = new_tracks(dd);
                    prev_track = terminated_tracks(dd);
                    if(tracks(prev_track, 3) == frameIdx-1)
                        continue;
                    end
                    tracks(prev_track, 3) = frameIdx-1;
                    tmp_L2(tmp_L2 == tracks(new_track, 1)) = tracks(prev_track, 1);
                    tracks(prev_track, 5) = 0;
                    tracks(prev_track, 6) = frameIdx-1;
                    tracks(new_track, :) = [-1, 0, 0, 0, 0, 0];
                end
            end
        end    
        
        tracks(tracks(:, 1) == -1, :) = [];
        % Labels (propagated)
        L2 = tmp_L2;
        L2_FOI(L2_FOI>0) = L2(L2_FOI>0);
    end
  end 
end