function [ Precision, Recall, F_Score ] = calcRecallPrecision( positives, TP, FP )
% Detection evaluation
% Signature:
%   [ Precision, Recall, F_Score ] = calcRecallPrecision( positives, TP, FP )
% Inputs:
%   positives: Number of positive ground truth elements.
%   TP: Number of true positive detections.
%   FP: Number of false positive detections.
% Outputs:
%   Precision, Recall, and F1 Score.
%%
%%
    FN = positives - TP;
    Precision = double(TP)/double(FP+TP);
    Recall = double(TP)/double(TP+FN);
    F_Score = 2*(Precision*Recall)/(Precision+Recall); 
end

