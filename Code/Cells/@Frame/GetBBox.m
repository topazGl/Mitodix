function bbox = GetBBox(image, COMyx)
% Get a bounding box for all non zero values in image.
%   bbox = GetBBox(image, COMyx)
% Inputs:
%   image = input image.
%   COMyx = (optional) point to center BBox around.
% Output:
%   bbox = struct:
%       bbox.y = [min row index, max row index];
%       bbox.x = [min col index, max col index];
    [y,x] = find(image > 0);
    if(isempty(y))
        error('image has no positive values');
    end
    
    bbox.y = [max(min(y)-1, 1), min(max(y)+1, size(image, 1))];
    bbox.x = [max(min(x)-1, 1), min(max(x)+1, size(image, 2))];
    
    % keep center of mass in the center of the box:
    if(exist('COMyx', 'var'))
        y = min(abs(COMyx(1) - bbox.y(1)), abs(bbox.y(2) - COMyx(1)));
        x = min(abs(COMyx(2) - bbox.x(1)), abs(bbox.x(2) - COMyx(2)));
        bbox.y = [COMyx(1)-y, COMyx(1)+y];
        bbox.x = [COMyx(2)-x, COMyx(2)+x];
    end
end

