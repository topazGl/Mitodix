function [allTPs, allFPs, TPR, FPR, precision] = ROC( mtdx, featureVec, featureName, allCandidates, GT, detectionRadii, autoThr, thrs, historyNum, dt )
    if(isempty(featureVec))
        allTPs = [];
        allFPs = [];
        TPR = nan;
        FPR = nan;
        precision = nan;
        return;
    end
    
    if(~exist('dt', 'var') || isempty(dt))
        dt = 1;
    end
    
    if(~exist('historyNum', 'var') || isempty(historyNum))
        historyNum = 2;
    end
    
    nCount = length(allCandidates);
    
    if(~exist('detectionRadii', 'var') || isempty(detectionRadii))
        lastFrame = mtdx.GetFrame(mtdx.LastFrameNum);
        detectionRadii = lastFrame.CellLength.Mean/2;
    end
    
    % output folder:
    outputFolder = fullfile(mtdx.OutputPath, 'patches');
    if(~exist(outputFolder, 'dir'))
        mkdir(outputFolder);
    end
            
    outputFolder = fullfile(outputFolder, featureName);
    if(~exist(outputFolder, 'dir'))
        mkdir(outputFolder);
    end
    
    outputPlots = fullfile(outputFolder, 'plots');
    if(~exist(outputPlots, 'dir'))
        mkdir(outputPlots);
    end
    
    if(~exist('GT', 'var') || isempty(GT))
        filtered_GT = mtdx.FilteredGT();
    else
        filtered_GT = GT;
    end
    GT_N = size(filtered_GT, 1);


    
    if(~exist('autoThr', 'var') || isempty(autoThr))
        % auto thr:
        [otsuThr, ~, ~, ~, ~, f] = AutoThr( featureVec, featureName, 1 );
        g = sum(featureVec < otsuThr)/nCount;
        if(g > 0.1)
            [otsuThr, ~, ~, ~, ~, f] = AutoThr( featureVec(featureVec < otsuThr), featureName );
        end
        try
            saveas(f, fullfile(outputPlots, [featureName, 'Sorted.jpg']));
            saveas(f, fullfile(outputPlots,[featureName, 'Sorted.fig']));
            close(f);
        end
        autoThr = otsuThr;
    end
    
    % divide into FP and TP:
    detections = checkDetections(filtered_GT, allCandidates, featureVec, detectionRadii, dt);
    allTPs = detections(detections ~= 0);
    featureTP = featureVec(allTPs, :);
    
    allFPs = ones(length(featureVec), 1);
    allFPs(allTPs) = 0;
    allFPs = find(allFPs == 1);
    featureFP = featureVec(allFPs, :);
    
    P = size(filtered_GT, 1);
    if(P > nCount)
        P = sum(detections ~= 0);
    end
    N = nCount - P;
    
    % feature histogram:
    [hFP, vFP] = hist(featureFP, 256);  % FP
    [hTP, vTP] = hist(featureTP, 256);  % TP

    % plot feature histogram:
    f = figure;
    bar(vFP, hFP); hold on;
    bar(vTP, hTP, 'FaceColor', 'r', 'EdgeColor', 'r', 'LineWidth', 0.9);
    dy = (max(hFP)-min(hFP))/length(vFP);
    ylim([min(hFP), max(hFP)]);
    ylabel('#Candidates');
    plot(autoThr*ones(length(vFP), 1), min(hFP):dy:dy*(length(vFP)-1), '--g', 'LineWidth', 3);
    hold off;
    l = legend({'FP', 'TP', 'auto thr'}, 'location', 'best');
    set(gca, 'FontSize', 40);
    set(l, 'FontSize', 40);
    saveas(f, fullfile(outputPlots, [featureName, 'Hist.jpg']));
    saveas(f, fullfile(outputPlots,[featureName, 'Hist.fig']));
    
    if(~isempty(featureTP))
        ylim([min(hTP), max(hTP)]);
        xlim([min(vTP), max(vTP)]);
    end
    
    saveas(f, fullfile(outputPlots, [featureName, 'HistZoom.jpg']));
    saveas(f, fullfile(outputPlots,[featureName, 'HistZoom.fig']));
    close(f);
    
    if(~GT_N)
        TPR = [];
        FPR = [];
        precision = [];
        return;
    end
    
    % thrs for ROC:
    if(~exist('thrs', 'var') || isempty(thrs))
        [~, v] = hist(featureVec, 256);
        thrs = v';
    end
    
    thrs = sort(unique([autoThr; thrs]));
    L = length(thrs);
    
    TPR = zeros(L, 1);
    FPR = zeros(L, 1);
    F_score = zeros(L, 1);
    precision =  zeros(L, 1);
    
    detectionCand = zeros(GT_N, 1);
    detected = zeros(GT_N, 1);
    
    for tt=1:L
        detectedAsT = find(featureVec < thrs(tt));
        T = length(detectedAsT);
        tmp_detections = checkDetections(filtered_GT, allCandidates(detectedAsT), featureVec(detectedAsT), detectionRadii, dt);
        TP = sum(tmp_detections~=0);
        FP = T - TP;
        TPR(tt) = TP/P;
        FPR(tt) = FP/N;
        [ precision(tt), ~, F_score(tt) ] = calcRecallPrecision( P, TP, FP );
        
        % Mark first detections:
        ind = find(~detected & tmp_detections~=0);
        for ii = 1:length(ind)
            i = ind(ii);
            detectionCand(i) = tmp_detections(i);
            detected(i) = thrs(tt);
        end
    end
    
    % detections:
    detected(detected == 0) = inf;
    [detected, ind] = sort(detected, 'descend'); % sort by detection thr
    detectionData.Cand = detectionCand(ind);
    detectionData.DetectionThr = detected;
    detectionData.GT = ind;
    
    idcs = detectionData.DetectionThr < autoThr + max(abs(featureVec))/10;
    highlight.Cand = detectionData.Cand(idcs);
    highlight.GT = ind(idcs);
    
%     if(exist('weights', 'var') && ~isempty(weights))
%         % plot "rejection" cells:
%         p = length(highlight.GT);
%         r = ceil(sqrt(p));
%         q = ceil(p/r);
%         f = figure;
%         for ii=1:p
%             if(highlight.Cand(ii) == 0)
%                 j = 0;
%             else
%                 j = highlight.Cand(ii);
%             end
%             if(j == 0)
%                 a = nan*ones(size(data, 2), 1);
%                 b = nan*ones(size(weights, 2), 1);
%             else
%                 a = data(j,:);
%                 b = weights(j,:);
%             end
%             b_im = reshape(b, [pSize, pSize]);
%             cell_im = [];
%             for hh=0:1:historyNum-1
%                 a_im = reshape(a(hh*sz+1:hh*sz + sz), [pSize, pSize]);
%                 cell_im = [cell_im; a_im];
%             end
%             cell_im = [cell_im; b_im/max(b)];
%             subplot(r,q,ii); title(['cell #',num2str(j)]); imagesc(cell_im);
%             if(ii == 1)
%                 colormap gray;
%             end
%             set(gca,'Yticklabel',[]); set(gca,'Xticklabel',[]);
%         end
%         saveas(f, fullfile(outputPlots, [featureName, 'Rejection.jpg']));
%         saveas(f, fullfile(outputPlots,[featureName, 'Rejection.fig']));
%         close(f);
%     end
    
    % plot ROC:
    [~, idx_maxRecall] = max(TPR);
    idx_maxRecall = min(idx_maxRecall);
    
    workPoint = autoThr;
    [~, idx_wp] = min(abs(thrs - autoThr));
    
    [~, idx_F] = max(F_score);
    idx_F = min(idx_F);
    
    f = figure;
    plot3(TPR, precision, thrs, '.-', 'linewidth', 1.5); xlabel('Recall (TPR)'); ylabel('Precision'); view(2); hold on;
    plot3(TPR(idx_wp), precision(idx_wp), workPoint, '*r', 'linewidth', 1.5);
    plot3(TPR(idx_F), precision(idx_F), thrs(idx_F), 'og', 'linewidth', 1.5); hold off;
    if(max(precision) && max(TPR))
        xlim([0,max(TPR)]);
        ylim([0,max(precision)]);
    end
    legend('Precision-Recall', ['auto calculated work point (R=', num2str(TPR(idx_wp), 3), ', P=',num2str(precision(idx_wp), 3),')'], ['max F-Score=', num2str(F_score(idx_F), 3), ' (R=', num2str(TPR(idx_F), 3), ', P=',num2str(precision(idx_F), 3),')'], 'Location', 'south');
    saveas(f, fullfile(outputPlots, [featureName, 'PRgraph.jpg']));
    saveas(f, fullfile(outputPlots,[featureName, 'PRgraph.fig']));
    close(f);
    
    f = figure;
    plot3(FPR, TPR, thrs, '.-', 'linewidth', 1.5); xlabel('FPR'); ylabel('Recall'); view(2); hold on;
    plot3(FPR(idx_wp), TPR(idx_wp), workPoint, '*r', 'linewidth', 1.5);
    plot3(FPR(idx_maxRecall), TPR(idx_maxRecall), thrs(idx_maxRecall), 'og', 'linewidth', 1.5); hold off;
    legend('ROC', ['auto calculated work point, R=', num2str(TPR(idx_wp), 3), ' (FPR=',num2str(FPR(idx_wp), 3),')'], ['max R=', num2str(TPR(idx_maxRecall), 3), ' (FPR=',num2str(FPR(idx_maxRecall), 3),')'], 'Location', 'south');
    xlim([0,0.5*min(1, max(FPR))]); ylim([0,1]);
    saveas(f, fullfile(outputPlots, [featureName, 'ROC.jpg']));
    saveas(f, fullfile(outputPlots,[featureName, 'ROC.fig']));
    close(f);
    
    save(fullfile(outputFolder,'analysis.mat'), 'autoThr', 'thrs', 'TPR', 'FPR', 'F_score', 'precision', 'detections', 'FP', 'TP', 'detectionData', 'highlight', 'allTPs', 'allFPs', '-v7.3');
    save(fullfile(outputFolder, 'filtered_GT.mat'), 'filtered_GT', '-v7.3');
    zip(fullfile(outputFolder,'Mitodix.zip'),fullfile('..\Mitodix'));
end

