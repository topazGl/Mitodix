function [newLabeled, split] = WaterShed( labeled, h, conn )
% Watershed for input labeled image.
%%
    if(~exist('conn', 'var') || isempty(conn))
        conn = 4;
    end
    
    if(~exist('h', 'var') || isempty(conn))
        h = 2;
    end
    
    bw = labeled~=0;
    D = -bwdist(~bw);
    mask = imextendedmin(D,h,conn);
    D2 = imimposemin(D,mask);
    newBW = watershed(D2);
    split = (newBW==0) & (bw);
%     newLabeled(~newBW) = 0;
%     newBW(~bw) = 0;
    
    newLabeled = labeled;
    newLabeled(split) = 0;
end

