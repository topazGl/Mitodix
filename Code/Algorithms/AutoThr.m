function [ otsuThr, h, v, sortedFeatureVec, idx, f ] = AutoThr( featureVec, featureName, plotFlag, C, outliersPer )
% Calculate threshold for the input feature vector and sort it.
% Signature:
%   [ otsuThr, h, v, sortedFeatureVec, idx, f ] = AutoThr( featureVec, featureName, plotFlag, C, outliersPer )
% Inputs:
%   featureVec: 1D feature vector.
%   featureName: used for plot title (relevant for plotFlag=1).
%   plotFlag: 0/1 flag (Default = 0).
%   C: C>1 = number of classes (Default = 20).
%   outliersPer: 2x1 vector of edge allowed percentage of outliers out of
%   the entire population (Default=[0.02, 0.08]).
% Outputs:
%   otsuThr: Otsu's threshold.
%   h: input feature histogram counts (256 bins).
%   v: input feature histogram bin levels (256 bins).
%   sortedFeatureVec: sorted feature vector.
%   idx: indices of input elements after sorting.
%   f: handle to output figure (relevant for plotFlag=1).
% Example:
%   I = imread('coins.png');
%   featureVec = I(:); % Intensity values as features.
%   [otsuThr, h, v, sortedFeatureVec, idx] = AutoThr( featureVec, 'Intensities', 1, 2 );
% See also: graythresh; sort
%%
    if(isempty(featureVec))
        otsuThr=1;
        h=nan;
        v=nan;
        sortedFeatureVec=[];
        idx=[];
        f=[];
        return;
    end
    
    % Initialize inputs:
    if(~exist('C', 'var') || isempty(C))
        C = 20;
    end
    
    if(~exist('plotFlag', 'var') || isempty(plotFlag))
        plotFlag = 0;
    end
    f = [];
    
    if(~exist('outliersPer', 'var') || isempty(outliersPer))
        outliersPer = [0.02, 0.08];
    end
    
    % Sort feature vector:
    [sortedFeatureVec, idx] = sort(double(featureVec));

    % Calculate thr:
    [h, v] = hist(double(featureVec), 256);
    otsuThr = v(1) + graythresh(featureVec)*(v(end)-v(1));
    
    ratio = sum(featureVec<otsuThr)/length(sortedFeatureVec);
    
    if(C>2)
        %4% outliers:
        if(ratio>outliersPer(end)/2)
            multi_thr = sort(multithresh(featureVec(featureVec<otsuThr),C-1));
            thresh = multi_thr(end);
            if(~isempty(thresh) && ~isempty(otsuThr) && ((sum(featureVec<thresh)/length(featureVec))>0.03))
                otsuThr = thresh;
            end
            ratio = sum(featureVec<otsuThr)/length(sortedFeatureVec);
            if((ratio>mean(outliersPer)) && (length(multi_thr)>1))
                otsuThr = multi_thr(end-1);
            end
        end
    end
    energy = cumsum(double(h)/sum(h));
    minV = v(find(energy >= outliersPer(1), 1, 'first'));
    maxV = v(find(energy <= outliersPer(end), 1, 'last'));
    otsuThr = min(max(otsuThr, minV), maxV);
    
    [~, dx_otsu] = min(abs(sortedFeatureVec - otsuThr));
    dx_otsu = dx_otsu(1);

    % plot:
    if(plotFlag)
        f = figure;
        subplot(2,1,1);
        bar(v, h);
        h_vec = min(h):max(h); hold on;
        plot(otsuThr*ones(size(h_vec)), h_vec, '--');
        hold off;
        title([featureName, ' histogram']);
        legend({['histogram ' featureName], ['auto thr=', num2str(otsuThr)]}, 'location', 'best');
        xlabel('bins values');
        ylabel('histogram');
        set(gca, 'FontSize', 20);
        
        smoothSorted = smooth(smooth(smooth(sortedFeatureVec)));
        subplot(2,1,2);
        plot(2:length(smoothSorted), smoothSorted(2:end), 'LineWidth', 2); hold on;
        plot(dx_otsu, otsuThr, '*g', 'LineWidth', 2);
        title(featureName);
        legend({['sorted ' featureName], 'auto thr'}, 'location', 'best');
        ylabel([featureName, ' value']);
        set(gca, 'FontSize', 20);
        ax1 = gca;
        
        set(f, 'Position', get(0, 'Screensize'));
    end
end

