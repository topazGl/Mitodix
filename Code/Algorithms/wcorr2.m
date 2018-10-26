function [wCorr, wMuRatio, wCov, wSim] = wcorr2( A, B, W, r, centerFlag )
% Weighted correlation and similarity.
% Signature:
%   [wCorr, wMuRatio, wCov, wSim] = wcorr2( A, B, W, r, centerFlag )
% Inputs:
%   A, B: are the two input matrices (same dimensions).
%   W: Weights matrix (same dimensions as A and B).
%   r: highest value percentage (Default = 1).
%   centerFlag: center A to mean 0 and B to mean 0 before calculation (Default = 1).
% Outputs:
%   wCorr: Weighted Pearson correlation coefficient.
%   wMuRatio: Ratio between weighted mean values.
%   wCov: Weighted covariance.
%   wSim: similarity.
% Example:
%     I = imread('coins.png');
%     A = I(75:130, 237:292); B = I(145:200,207:262);
%     BW_A = imfill(imclearborder(im2bw(A)),'holes');
%     BW_B = imfill(imclearborder(im2bw(B)),'holes');
%     BW = BW_A | BW_B;
%     W = CalcPHI(CalcSDF(BW));
%     [wCorr, wMuRatio, wCov, wSim] = wcorr2( A, B, W, [], 0);
%     disp(['wCorr = ', num2str(wCorr), '; wSim = ', num2str(wSim), '; corr2 = ', num2str(corr2(A, B))]);
%     figure; 
%     subplot(3,1,1); imagesc(A); title('Image A'); colorbar; ax1 = gca; cc = get(gca, 'clim');
%     subplot(3,1,2); imagesc(B); title('Image B'); colorbar; ax2 = gca; set(gca, 'clim', cc);
%     subplot(3,1,3); imagesc(W/sum(W(:))); title('Weights'); colormap jet;
%     colorbar; linkaxes([ax1, ax2, gca]);
%
% See also: CalcPHI, CalcSDF, corr2
%
% Author: T. Gilad, 2018
%%
    if(~isequal(size(A), size(B)))
        error('A and B must be of same size!');
    end
    
    if(~exist('W', 'var') || isempty(W) || isequal(W, 1))
        W = ones(size(A));
        compVal = corr2(A, B);
    else
        compVal = [];
        if(~isequal(size(A), size(W)))
            error('A and W must be of same size!');
        end
    end
    
    if(~exist('centerFlag', 'var') || isempty(centerFlag))
        centerFlag = 1;
    end
    
    if(~exist('r', 'var') || isempty(r))
        r = 1;
    else
        r = max(0, r);
        r = min(1, r);
    end
    
    % double precision:
    A = double(A);
    B = double(B);
    W = double(W);
    
%     figure; subplot(3,1,1); imagesc(A.*W); cc = get(gca, 'clim'); colorbar; subplot(3,1,2); imagesc(B.*W); set(gca, 'clim', cc); colorbar; subplot(3,1,3); imagesc(W/max(W(:))); colorbar;
    
    % weights must sum to 1:
    W = W / sum(sum(W));
    
    % vector representation:
    a = A(:);
    b = B(:);
    w = W(:);
    
    % weighted mean:
    wMuA = a'*w;
    wMuB = b'*w;
    
    wMuRatio = muRatio(a, b, w, r, wMuA, wMuB);
    [~, histSim] = weightedHist(a, b, w);
    
    % center data:
    if(centerFlag)
        a = a - wMuA;
        b = b - wMuB;
    end
    
    % weighted covariance:
    wAutoCovA = a'*(w.*a);
    wAutoCovB = b'*(w.*b);
    wCov = a'*(w.*b);
    
    % weighted correlation coeff:
    wCorr = wCov/sqrt(wAutoCovA * wAutoCovB);
    
    wSim =  0.5*(wCorr+1) * wMuRatio ;

    % competability to corr2 in case of no weights:
    if(~isempty(compVal))
        compDiff = abs(compVal - wCorr);
        if(compDiff > 10^-5)
            warning(['corr2 and wcorr2 with W=1 diff = ', num2str(compDiff)]);
        end
    end
end

function wMuRatio = muRatio(a, b, w, r, wMuA, wMuB)
    % weighted mean ratio:
    if(r < 1 && r > 0)
        % Consider only highest intensity values:
        wA = w.*a;
        [~, v] = hist(wA, 256);
        p = floor((1-r)*length(v));
        ind = wA > v(p);
        rA = (a(ind)'*w(ind))/sum(w(ind));
        
        wB = w.*b;
        [~, v] = hist(wB, 256);
        p = floor((1-r)*length(v));
        ind = wB > v(p);
        rB = (b(ind)'*w(ind))/sum(w(ind));
    else
        rA = wMuA;
        rB = wMuB;
    end
    
    if(rB == 0)
        wMuRatio = (rA == 0);
    else
        wMuRatio = rA / rB;
    end
    if(wMuRatio > 1)
        wMuRatio = 1/wMuRatio;
    end
end

function [muRatio, histSim, bins, ha, hb] = weightedHist(a, b, w, nbins)
    if(~exist('nbins', 'var') || isempty(nbins) || nbins < 1)
        nbins = 50;
    elseif(nbins < 5)
        nbins = 5;
    end
    w = w/max(w);
    [~, va] = hist(a.*w, nbins);
    [~, vb] = hist(b.*w, nbins);
    v_min = min([va, vb]);
    v_max = max([va, vb]);
    binWidth = (v_max-v_min)/nbins;
    v = v_min:binWidth:v_max;%unique(sort([va, vb]))';
    ha = zeros(length(v)-1, 1);
    hb = zeros(length(v)-1, 1);
    bins = zeros(length(v)-1, 1);
    w = w/sum(w);
    for ii=2:length(v)
        aIndx = (a > v(ii-1) & (a <= v(ii)));
        ha(ii-1) = a(aIndx)'*w(aIndx);
        bIndx = (b > v(ii-1) & (b <= v(ii)));
        hb(ii-1) = b(bIndx)'*w(bIndx);
        bins(ii-1) = 0.5*(v(ii-1) + v(ii));
    end
    histSim = 2*sqrt(ha'*hb)/(sqrt(ha'*ha)+sqrt(hb'*hb));
    muRatio = ((bins'*ha)/sum(ha))/((bins'*hb)/sum(hb));
    if(muRatio > 1)
        muRatio = 1 / muRatio;
    end
%     f = figure; subplot(2,1,1); bar(bins, ha); subplot(2,1,2); bar(bins, hb); [histSim, muRatio]
%     close(f);
end