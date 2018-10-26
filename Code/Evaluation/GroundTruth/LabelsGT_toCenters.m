function [ GT ] = LabelsGT_toCenters( GT, mtdx )
    N = size(GT, 1);
    GT = table2struct(GT);
    for n=1:N
        motherFrame = mtdx.GetFrame(GT(n).BirthFrame-1);
        if(~isempty(motherFrame))
            c = motherFrame.GetCell(GT(n).MotherLabel);
            if(~isempty(c))
                GT(n).MotherCOM = c.COM;
            end
        end
        
        birthFrame = mtdx.GetFrame(GT(n).BirthFrame);
        if(~isempty(birthFrame))
            c = birthFrame.GetCell(GT(n).Daughter1Label);
            if(~isempty(c))
                GT(n).Daughter1COM = c.COM;
            end
            
            c = birthFrame.GetCell(GT(n).Daughter2Label);
            if(~isempty(c))
                GT(n).Daughter2COM = c.COM;
            end
        end
    end
    GT = struct2table(GT);
end

