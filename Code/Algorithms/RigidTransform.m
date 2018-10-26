function im = RigidTransform(im, dt, d_deg, isMask)
    if(~exist('isMask', 'var') || isempty(isMask))
        isMask = 0;
    end
    
    se = strel('disk',1);
    RA = imref2d(size(im));
    if(dt)
        im = imtranslate(im, RA, dt, 'cubic');
        if(isMask)
            im = imclose((im ~= 0), se);
        end
    end
    
    if(d_deg)
        im = imrotate(im, d_deg, 'bicubic', 'crop');
        if(isMask)
            im = imclose((im ~= 0), se);
        end
    end
end