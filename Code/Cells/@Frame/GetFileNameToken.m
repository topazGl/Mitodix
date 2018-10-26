function [t, fileNames, s] = GetFileNameToken( fileNames, expr )
% Get file names token according to input expression.
%   [t, fileNames] = GetFileNameToken( fileNames, expr )
% Inputs:
%   fileNames = cell array of file names.
%   expr = token expression in file name.
% Outputs:
%   t = file names tokens.
%   fileNames = all file names with a token.
%     if(ischar(expr))
%         expr = {expr};
%     end
    if(ischar(fileNames))
        fileNames = {fileNames};
    end
    
    tokens = cellfun(@(str) regexp(str,expr,'tokens'),fileNames,'uniformoutput',false);
    validFiles = cellfun(@(t)~isempty(t),tokens);
    tokens = tokens(validFiles);
    fileNames = fileNames(validFiles);
    s = cellfun(@(token) str2double(token{1}{1}),tokens);
    t = cellfun(@(token) str2double(token{1}{end}),tokens);
    [~, ind] = sort(t);
    fileNames = fileNames(ind);
    s = s(ind);
    [~, ind] = sort(s);
    fileNames = fileNames(ind);
    t = t(ind);
end

