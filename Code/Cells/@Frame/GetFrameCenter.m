function frameCenterPoint = GetFrameCenter(X, Y)
% Calculate the middle pixel of a frame the size [X, Y].
%   frameCenterPoint = GetFrameCenter(X, Y)
% Inputs:
%   X = frame X size [pixels].
%   Y = frame Y size [pixels].
% Output:
%   frameCenterPoint = frame [x,y] center point.

    frameCenterPoint = [0, 0];

    if(mod(X, 2) == 0)
        % even number of columns:
        frameCenterPoint(1) = X/2;
    else
        % odd number of columns:
        frameCenterPoint(1) = (X-1)/2;
    end

    if(mod(Y, 2) == 0)
        % even number of rows:
        frameCenterPoint(2) = Y/2;
    else
        % odd number of rows:
        frameCenterPoint(2) = (Y-1)/2;
    end
end

