function [ performance, falsePositiveIdxs, falseNegativeIdxs ] = evaluation( mtdx, filtered_GT, candidates, feature, allN, radi, timeDiff )
    if(~exist('radi', 'var') || isempty(radi))
        lastFrame = mtdx.LastFrame;
        radi = lastFrame.CellLength.Mean/2;
    end

    if(~exist('timeDiff', 'var') || isempty(timeDiff))
        timeDiff = 1;
    end
    
    performance = cell(3,1);
    
    P = size(filtered_GT, 1);
    N = allN - P;
    T = size(candidates, 1);
    if(T==1)
        T = length(candidates);
    end
    
    i=1;
    detections = checkDetections(filtered_GT, candidates, feature, radi, timeDiff, 'all');
    
%     FNidx = find(detections <= 0);
%     for kk=1:length(FNidx)
%         fNum = filtered_GT(FNidx(kk),:).BirthFrame;
%         f = figure; birthFrame = mtdx.GetFrame(fNum); motherFrame = mtdx.GetFrame(fNum-1); subplot(2,2,1); imagesc(motherFrame.Image); cc = get(gca, 'clim'); ax1 = gca; subplot(2,2,2); imagesc(motherFrame.Labeled); ax2 = gca; subplot(2,2,3); imagesc(birthFrame.Image); set(gca, 'clim', cc); ax3 = gca; subplot(2,2,4); imagesc(birthFrame.Labeled); linkaxes([ax1, ax2, ax3, gca]);
%         close(f);
%     end

    TP = sum(detections>0);
    FP = T - TP;
    performance{i}.TP = TP;
    performance{i}.FP = FP;
    performance{i}.FN = P-TP;
    performance{i}.TPR = TP/P;
    performance{i}.precision = TP/T;
    performance{i}.FPR = FP/allN;
    performance{i}.Fscore = 2*(performance{i}.precision*performance{i}.TPR)/(performance{i}.precision+performance{i}.TPR);
    
    falsePositiveIdxs = [];
    for cc=1:length(candidates)
        if(~any(detections == cc))
            falsePositiveIdxs = [falsePositiveIdxs; cc];
        end
    end
    falseNegativeIdxs = find(detections==0);
    
    i = i+1;
    detections = checkDetections(filtered_GT, candidates, feature, radi, timeDiff, 'mothers');
    TP = sum(detections>0);
    FP = T - TP;
    performance{i}.TP = TP;
    performance{i}.FP = FP;
    performance{i}.FN = P-TP;
    performance{i}.TPR = TP/P;
    performance{i}.precision = TP/T;
    performance{i}.FPR = FP/allN;
    performance{i}.Fscore = 2*(performance{i}.precision*performance{i}.TPR)/(performance{i}.precision+performance{i}.TPR);
    
    i = i+1;
    detections = checkDetections(filtered_GT, candidates, feature, radi, timeDiff, 'daughters');
    TP = sum(detections>0);
    FP = T - TP;
    performance{i}.TP = TP;
    performance{i}.FP = FP;
    performance{i}.FN = P-TP;
    performance{i}.TPR = TP/P;
    performance{i}.precision = TP/T;
    performance{i}.FPR = FP/allN;
    performance{i}.Fscore = 2*(performance{i}.precision*performance{i}.TPR)/(performance{i}.precision+performance{i}.TPR);
end

