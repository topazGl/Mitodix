function [detections, connections, TPs] = checkDetections(filtered_GT, candidates, feature, radi, timeDiff, checkType )
    N = size(filtered_GT, 1);
    if(~N)
        detections = [];
        connections = [];
        TPs = [];
        return;
    end
    if(iscell(candidates))
        candidates = cell2mat(candidates);
    end
    
    daughter1_only = 0;
    detections = zeros(N,1);
    connections = cell(N,1);
    
    if(isempty(candidates))
        return;
    end
    
    if(~exist('checkType', 'var') || isempty(checkType))
        checkType = 'all';
    else
        checkType = lower(checkType);
    end
    
    if(~exist('radi', 'var') || isempty(radi))
        radi = 20;
    end
    
    if(~exist('timeDiff', 'var') || isempty(timeDiff))
        timeDiff = 1;
    end

    if(~exist('feature', 'var') || isempty(feature))
        if(isempty(candidates))
            return;
        elseif (isfield(candidates, 'RegMatch'))
            feature = double(cat(1, candidates.RegMatch));
        else
            feature = ones(size(candidates,1), 1);
        end
    end

    if(~isequal(checkType, 'daughters'))
        if(isfield(candidates, 'COM'))
        mothersCenters = cat(1, candidates.COM);
        mothersFlag = 1;
        elseif(isfield(candidates, 'MotherCOM'))
            mothersCenters = cat(1, candidates.MotherCOM);
            mothersFlag = 1;
            if(any(mothersCenters(:) == -1))
                mothersFlag = 0;
                checkType = 'daughters';
            end
        else
            mothersFlag = 0;
        end
    else
        mothersFlag = 0;
    end
    
    if(~isequal(checkType, 'mothers') && isfield(candidates, 'Daughter1COM'))
        centers1 = cat(1, candidates.Daughter1COM);
        centers2 = cat(1, candidates.Daughter2COM);
        daughter1_only = any(centers2(:) == -1);
        daughtersFlag = 1;
    else
        daughtersFlag = 0;
    end
    
    if(isfield(candidates, 'FrameNum'))
        frames = double(cat(1, candidates.FrameNum));
    else
        frames = double(cat(1, candidates.BirthFrame));
    end
    
    for n=1:N
        timeDist = abs(frames - double(filtered_GT(n,:).BirthFrame));
        candIdx = find(timeDist <= timeDiff);
        if(isempty(candIdx))
            continue;
        end
        
        if(mothersFlag)
            subCenters = mothersCenters(candIdx, :);
%             if(isequal(checkType, 'one_daughter'))
%                 dists1 = sqrt((subCenters(:,1) - filtered_GT(n,:).Daughter1COM(1)).^2 + (subCenters(:,2) - filtered_GT(n,:).Daughter1COM(2)).^2);
%                 dists2 = sqrt((subCenters(:,1) - filtered_GT(n,:).Daughter2COM(1)).^2 + (subCenters(:,2) - filtered_GT(n,:).Daughter2COM(2)).^2);
%                 dists = min([dists1, dists2]);
%             else
                dists = sqrt((subCenters(:,1) - filtered_GT(n,:).MotherCOM(1)).^2 + (subCenters(:,2) - filtered_GT(n,:).MotherCOM(2)).^2);
%             end
            motherCondition = dists <= radi;
        else
            motherCondition = ones(length(candIdx), 1);
        end
        
        if(daughtersFlag)
            subCenters1 = centers1(candIdx, :);
            subCenters2 = centers2(candIdx, :);
            
            dists11 = sqrt((subCenters1(:,1) - filtered_GT(n,:).Daughter1COM(1)).^2 + (subCenters1(:,2) - filtered_GT(n,:).Daughter1COM(2)).^2);
            dists22 = sqrt((subCenters2(:,1) - filtered_GT(n,:).Daughter2COM(1)).^2 + (subCenters2(:,2) - filtered_GT(n,:).Daughter2COM(2)).^2);
            
            dists12 = sqrt((subCenters1(:,1) - filtered_GT(n,:).Daughter2COM(1)).^2 + (subCenters1(:,2) - filtered_GT(n,:).Daughter2COM(2)).^2);
            dists21 = sqrt((subCenters2(:,1) - filtered_GT(n,:).Daughter1COM(1)).^2 + (subCenters2(:,2) - filtered_GT(n,:).Daughter1COM(2)).^2);
            
            % One of the cells was not annotated:
            if(filtered_GT(n,:).Daughter1COM(1) < 0)
                dists11 = zeros(length(dists11), 1);
                dists21 = zeros(length(dists21), 1);
            end
            if(filtered_GT(n,:).Daughter2COM(1) < 0)
                dists12 = zeros(length(dists12), 1);
                dists22 = zeros(length(dists22), 1);
            end
            if(daughter1_only)
                dists21 = zeros(length(dists21), 1);
                dists22 = zeros(length(dists22), 1);
            end
            daughtersCondition = ((dists11 <= radi) & (dists22 <= radi)) | ((dists12 <= radi) & (dists21 <= radi));
        else
            daughtersCondition = ones(length(candIdx), 1);
        end
        
        candIdx = candIdx(motherCondition & daughtersCondition);
        if(isempty(candIdx))
            continue;
        end
        
        connections{n} = candIdx;
        
%         connections{n} = find(dists<=radi);
%         fNums{n} = double(frames(connections{n}));
%         idx = find(abs(fNums{n} - filtered_GT(n,:).BirthFrame) <= timeDiff);
        
        % Each candidate can only be connected to one GT:
        if(~isempty(candIdx) && any(detections~=0))
            detectedCands = detections(detections~=0);
            filteredIdx = [];
            for jj = 1:length(candIdx)
                if ~any(detectedCands == candIdx(jj))
                    filteredIdx =  [filteredIdx; candIdx(jj)];
                end
            end
            candIdx = filteredIdx;
        end
        
        foundFlag = ~isempty(candIdx);
        detections(n) = foundFlag;
        if(~foundFlag)
            continue;
        end
        
%         % Match is the closest candidate COM inside the smallest time range:
%         timeDist = timeDist(candIdx);
%         candIdx = candIdx(timeDist == min(timeDist));
%         
%         dists = dists(candIdx);
%         [~, indcs] = sort(dists);

        % Match is the smallest cost inside the smallest time range:
        cost = feature(candIdx);
        candIdx = candIdx(cost == min(cost));
        
        timeDist = timeDist(candIdx);
        [~, indcs] = sort(timeDist);
        
        candIdx = candIdx(indcs);
        
        detections(n) = candIdx(1);
    end
    
    TPs = detections(detections ~= 0);
    if(~isempty(TPs) && length(TPs) > length(unique(TPs)))
        error('Some of the candidates are connected to more than one GT!');
    end
end