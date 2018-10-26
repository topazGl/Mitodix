function [index, found] = FindDoubleSortedElementPosition(array, val)
% Find an element in a double sorted array.
%   [index, found] = FindDoubleSortedElementPosition(array, val)
% Inputs:
%   array = input array.
%   val = value to search for.
% Outputs:
%   index = val index inside the array (0 if not found).
%   found = 1 if val was found, else 0.
    index = 0;
    found = 0;
    N = size(array, 1);
    while ((index < N) && (array(index+1,1) < val(1)))
        index = index+1;
    end
    while ((index < N) && (array(index+1,1) == val(1)) && (array(index+1,2) < val(2)))
      index = index+1;
    end
    index = index+1;
    if(index > 0 && index <= N)
        found = isequal(array(index,:), val);
    end
end

