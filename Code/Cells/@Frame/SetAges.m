function SetAges(self, ages)
% Set cells ages.
%   SetAges(self, ages)
% Inputs:
%   self = Frame object.
%   ages = array of input ages.
    try
        if (length(ages) ~= length(self.labels))
            throw('ages input vector must be the same size as labeles vector.');
        else
            self.ages = ages;
        end
    catch EX
        self.Log.mlog = {LogType.ERROR ...
                      ,mfilename('class') ...
                      ,[self.Log.Me, EX.message]};
        rethrow(EX);
    end
end

