function [performance] = evalTRA(path_TRA, path_GT, timeDiff, radi, allN)
% Evaluate tracking mitosis detection compared to ground truth
    if(~exist('radi', 'var') || isempty(radi))
        radi = 20;
    end

    if(~exist('timeDiff', 'var') || isempty(timeDiff))
        timeDiff = 1;
    end

    if(~exist('allN', 'var') || isempty(allN))
        allN = nan;
    end
    
    if(~exist(fullfile(path_TRA, 'mitosisTRA.mat'), 'file'))
        [ mitosisTRA ] = ReadISBI_GT(path_TRA, path_TRA);
        mitosisTRA = convertISBIGTformat( mitosisTRA, path_TRA );
        save(fullfile(path_TRA, 'mitosisTRA.mat'), 'mitosisTRA');
    else
        load(fullfile(path_TRA, 'mitosisTRA.mat'));
    end

    load(fullfile(path_GT, 'mitosisGT.mat'))

    detections = checkDetections(GT, mitosisTRA, [], radi, timeDiff );
    P = size(GT, 1);
    T = size(mitosisTRA, 1);
    performance.TP = sum(detections>0);
    performance.FP = T - performance.TP;
    performance.FN = P-performance.TP;
    performance.TPR = performance.TP/P;
    performance.precision = performance.TP/T;
    performance.FPR = performance.FP/allN;
    performance.Fscore = 2*(performance.precision*performance.TPR)/(performance.precision+performance.TPR);
end