function [detection_labels, connections] = labelTriplets(filtered_GT, candidates, radi, timeDiff )
    % 0 = all cells in triplets are false
    % 1 = 1 daughter + mother are false (one true daughter)
    % 2 = both daughters are false (true mother)
    % 3 = false daughter (both mother and one daughter are true)
    % 4 = false mother (both daughters are true)
    % 5 = all cells in triplet are true
    
    if(~exist('radi', 'var') || isempty(radi))
        radi = 20;
    end
    
    if(~exist('timeDiff', 'var') || isempty(timeDiff))
        timeDiff = 1;
    end
    
    N = size(filtered_GT, 1);
    
    if(isempty(candidates))
        connections = [];
        detection_labels = [];
        return;
    end
    
    M = length(candidates);
    detection_labels = zeros(M, 1);
    connections = zeros(M, 1);
    
    mothersCenters = cat(1, candidates.COM);
    centers1 = cat(1, candidates.Daughter1COM);
    centers2 = cat(1, candidates.Daughter2COM);
    frames = double(cat(1, candidates.FrameNum));

    for n=1:N
        timeDist = abs(frames - double(filtered_GT(n,:).BirthFrame));
        candIdx = find(timeDist <= timeDiff);
        if(isempty(candIdx))
            continue;
        end
        
        subCenters = mothersCenters(candIdx, :);
        dists = sqrt((subCenters(:,1) - filtered_GT(n,:).MotherCOM(1)).^2 + (subCenters(:,2) - filtered_GT(n,:).MotherCOM(2)).^2);
        motherCondition = dists <= radi;

        subCenters1 = centers1(candIdx, :);
        subCenters2 = centers2(candIdx, :);
        
        if(filtered_GT(n,:).Daughter1COM(1) < 0)
            dists11 = zeros(length(candIdx), 1);
        else
            dists11 = sqrt((subCenters1(:,1) - filtered_GT(n,:).Daughter1COM(1)).^2 + (subCenters1(:,2) - filtered_GT(n,:).Daughter1COM(2)).^2);
        end

        if(filtered_GT(n,:).Daughter2COM(1) < 0)
            dists12 = zeros(length(candIdx), 1);
        else
            dists12 = sqrt((subCenters1(:,1) - filtered_GT(n,:).Daughter2COM(1)).^2 + (subCenters1(:,2) - filtered_GT(n,:).Daughter2COM(2)).^2);
        end

        if(filtered_GT(n,:).Daughter1COM(1) < 0)
            dists21 = zeros(length(candIdx), 1);
        else
            dists21 = sqrt((subCenters2(:,1) - filtered_GT(n,:).Daughter1COM(1)).^2 + (subCenters2(:,2) - filtered_GT(n,:).Daughter1COM(2)).^2);
        end

        if(filtered_GT(n,:).Daughter2COM(1) < 0)
            dists22 = zeros(length(candIdx), 1);
        else
            dists22 = sqrt((subCenters2(:,1) - filtered_GT(n,:).Daughter2COM(1)).^2 + (subCenters2(:,2) - filtered_GT(n,:).Daughter2COM(2)).^2);
        end

        daughters_cond = (dists11 <= radi & dists22 <= radi) | (dists12 <= radi & dists21 <= radi);
        daughter1_cond = (dists11 <= radi) | (dists21 <= radi);
        daughter2_cond = (dists22 <= radi) | (dists12 <= radi);

        daughter1_cond_only =  daughter1_cond & (~daughter2_cond);
        daughter2_cond_only =  daughter2_cond & (~daughter1_cond);
        one_daughter_cond = xor(daughter1_cond_only, daughter2_cond_only);
        
        all_cond = (motherCondition & daughters_cond);
        daughters_only_cond = (daughters_cond & (~motherCondition));
        one_daughter_only_cond = one_daughter_cond & (~motherCondition);
        mother_and_one_daughter_cond = (motherCondition & one_daughter_cond);
        mother_only_cond = (motherCondition & (~mother_and_one_daughter_cond) & (~all_cond));

        % only true daughter (false mother and daugher):
        cond = one_daughter_only_cond;
        idx = candIdx(cond);
        [detection_labels(idx), connections(idx)] = set_label(idx, 1, n, detection_labels(idx), connections(idx));

                
        % only true mother (false daughters):
        cond = mother_only_cond;
        idx = candIdx(cond);
        [detection_labels(idx), connections(idx)] = set_label(idx, 2, n, detection_labels(idx), connections(idx));
        
        
        % true mother + daughter (one daughter false):
        cond = mother_and_one_daughter_cond;
        idx = candIdx(cond);
        [detection_labels(idx), connections(idx)] = set_label(idx, 3, n, detection_labels(idx), connections(idx));

        
        % only true daughters (false mother):
        cond = daughters_only_cond;
        idx = candIdx(cond);
        [detection_labels(idx), connections(idx)] = set_label(idx, 4, n, detection_labels(idx), connections(idx));
        
        % true mother and daughters:
        cond = all_cond;
        detection_labels(candIdx(cond)) = 5;
        connections(candIdx(cond)) = n;
    end
end

function [new_labels, new_connections] = set_label(idx, new_label, new_conn, current_labels, current_connections)  
    new_labels = current_labels;
    new_connections = current_connections;
    
    K = length(idx);
    for kk=1:K
        if(current_labels(kk) < new_label)
            new_labels(kk) = new_label;
            new_connections(kk) = new_conn;
        end
    end
end