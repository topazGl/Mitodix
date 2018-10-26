function orientation = CalcOrientation(dA, dB)
% Angle between the line that connects the two COMs and the frame's Y axis:
% Assumption: daughter1 has smaller Y coordinates than daughter2.
    dy = dA(2)-dB(2);
    dx = dA(1)-dB(1);

    if(dy == 0)
        orientation = -sign(dx)*90;
    elseif(dx == 0)
        orientation = (dy>0)*180;
    else
        orientation = sign(dx)*(90+sign(dy)*abs(atand(dy/dx)));
%                 if(dx < 0)
%                     orientation = orientation + sign(orientation)*180;
%                 end
    end

    if(isnan(orientation))
        orientation = 0;
        return;
    end

%             if(orientation < 0)
%                 orientation = 360 + orientation;
%             end

%             orientation = mod(orientation, 360);           
end