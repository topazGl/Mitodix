function RescaleGrayLevel(self, factor, uintFlag, offset)
% Rescale image Gray level to uint16
%   RescaleGrayLevel(self, factor)
% Inputs:
%   self = Frame object.
%   factor = intensity factor.
    if(~exist('uintFlag', 'var') || isempty(uintFlag))
        uintFlag = 1;
    end
    if(~exist('offset', 'var') || isempty(offset))
        offset = 0;
    end
    self.Image = (double(self.Image)+offset)*factor;
    self.Image(self.Image<0) = 0;
    
    if(uintFlag)
        self.Image = uint8(self.Image);
    end
end

