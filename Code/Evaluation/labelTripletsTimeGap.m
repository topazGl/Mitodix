function [detection_labels, connections] = labelTripletsTimeGap(filtered_GT, candidates, radi )
    % -1 = false
    % 0 = true in dt=0
    % 1 = true in dt=1
    % 2 = true in dt=2
    % 3 = true in dt=3
    % 4 = true in dt=4
    
    if(~exist('radi', 'var') || isempty(radi))
        radi = 20;
    end
    
    N = size(filtered_GT, 1);
    
    if(isempty(candidates))
        connections = [];
        detection_labels = [];
        return;
    end
    
    M = length(candidates);
    detection_labels = -1*ones(M, 1);
    connections = zeros(M, 1);
    
    mothersCenters = cat(1, candidates.COM);
    centers1 = cat(1, candidates.Daughter1COM);
    centers2 = cat(1, candidates.Daughter2COM);
    frames = double(cat(1, candidates.FrameNum));

    for timeDiff=4:-1:0
        for n=1:N
            timeDist = abs(frames - double(filtered_GT(n,:).BirthFrame));
            candIdx = find(timeDist <= timeDiff);
            if(isempty(candIdx))
                continue;
            end

            subCenters = mothersCenters(candIdx, :);
            dists = sqrt((subCenters(:,1) - filtered_GT(n,:).MotherCOM(1)).^2 + (subCenters(:,2) - filtered_GT(n,:).MotherCOM(2)).^2);
            motherCondition = dists <= radi + (5*timeDiff);

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

            % true mother and daughters:
            cond = (motherCondition & daughters_cond);
            idx = candIdx(cond);
            [detection_labels(idx), connections(idx)] = set_label(idx, timeDiff, n, detection_labels(idx), connections(idx));
        end
    end
end

function [new_labels, new_connections] = set_label(idx, new_label, new_conn, current_labels, current_connections)  
    new_labels = current_labels;
    new_connections = current_connections;
    
    K = length(idx);
    for kk=1:K
        if((current_labels(kk) < 0) || (current_labels(kk) > new_label))
            new_labels(kk) = new_label;
            new_connections(kk) = new_conn;
        end
    end
end
