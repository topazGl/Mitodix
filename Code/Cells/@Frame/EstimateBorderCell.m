function [isOnBorder, img, lbld] = EstimateBorderCell(self, cellLabel)
% Since cell is on the border of the frame,
% estimate its missing part assuming symmetry along
% major axis.
%   [isOnBorder, img, lbld] = EstimateBorderCell(self, cellLabel)
% Inputs:
%   self = Frame object.
%   cellLabel = cell label.
% Outputs:
%   isOnBorder = cell is on the border of the image flag.
%   img = farme gray-scale image with estimated cell.
%   lbld = farme labeled image with estimated cell.
    isOnBorder = self.IsLabelOnBorder(cellLabel);
    if(isOnBorder)
        org_img = self.PaddedImage;
        org_msk = Frame.CalcMask(self.PaddedLabeled == cellLabel);
        if(~any(org_msk(:)))
            error(['Label ', num2str(cellLabel), ' does not exist on frame ', num2str(self.Number), '.']);
        end
        
        % Since cell is on the border of the frame,
        % estimate its missing part assumin symmetry along
        % major axis:
        s = regionprops(org_msk,'Centroid');
        com = cat(1, s(1).Centroid);

        translation = self.Center - com;
        s = regionprops(org_msk,'Orientation');
        ang = (-1)*cat(1, s(1).Orientation);
        img = imtranslate(org_img,translation);
        img = imrotate(img, ang, 'bilinear', 'crop');
        msk = imtranslate(org_msk,translation);
        msk = imrotate(msk, ang, 'bilinear', 'crop');
        s = regionprops(msk,'Centroid');
        center = cat(1, s(1).Centroid);

        % mirror left-right:
        msk_flip = flipud(msk);
        img_flip = flipud(img);
        s = regionprops(msk_flip,'Centroid');
        center_flip = cat(1, s(1).Centroid);

        % center mirrored image:
        translation = center - center_flip;
        img_flip = imtranslate(img_flip,translation);
        msk_flip = imtranslate(msk_flip,translation);

        % Change back to original orientation:
        img_flip = imrotate(img_flip, -ang, 'bilinear', 'crop');
        msk_flip = imrotate(msk_flip, -ang, 'bilinear', 'crop');

        s = regionprops(msk_flip,'Centroid');
        center_flip = cat(1, s(1).Centroid);
        translation = com - center_flip;
        img_flip = imtranslate(img_flip,translation);
        msk_flip = imtranslate(msk_flip,translation);

        X = self.Size(1); Y = self.Size(2);
        paddedMargin = ones(X + 2*self.W, Y + 2*self.W);
        paddedMargin(self.W+1:self.W+X, self.W+1:self.W+Y) = 0;
            
        missing_msk = (msk_flip & not(org_msk) & paddedMargin);
        msk = (org_msk | missing_msk);
        img = org_img;
        img(missing_msk) = img_flip(missing_msk);

        % Make sure the result is only one blob:
        parts = double(labelmatrix(bwconncomp(msk)));
        if(length(parts) > 1)
            s = regionprops(parts,'Area');
            size = cat(1, s.Area);
            id = find(size == max(size), 1);
            msk = (parts == id);
        end

        lbld = self.PaddedLabeled;
        lbld(self.PaddedLabeled == cellLabel) = 0;
        lbld(msk) = cellLabel;

    else
        lbld = self.PaddedLabeled;
        img = self.PaddedImage;
    end           
end

