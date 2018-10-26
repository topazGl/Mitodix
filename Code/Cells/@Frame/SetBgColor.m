function SetBgColor(self, bgColor)
% Set back ground gray level color of frame image.
%   SetBgColor(self, bgColor)
% Inputs:
%   self = Frame object.
%   bgColor = new back ground color.
    self.Image = uint8(double(self.Image) - self.BgColor + bgColor);
    self.PaddedImage = uint8(double(self.PaddedImage) - self.BgColor + bgColor);
    self.BgColor = bgColor;
end

