classdef Statistics
    % statistical properties
   
    properties(Constant)
        Log = SingleInstance.Logger;     % Application logger.
    end
    
    properties (GetAccess = private)
        minV = 0;
        maxV = 0;
        meanV = 0;
        stdV = 0;
    end
    
    properties
        Min;
        Max;
        Mean;
        STD;
    end
    
    methods
        %% Constructor:
        function self = Statistics(values)
        % values is a vector of values.
            [m,n] = size(values);
            l = length(values(:));
            % verify that the input is a vector (not a matrix):
            if ((l > m) && (l > n))
                self.Log.mlog = {LogType.ERROR ...
                                  ,mfilename('class') ...
                                  ,[self.Log.Me,' Input values is a matrix instead of a vector.']};
                throw('In Statistics constructor: Input values is a matrix instead of a vector.');
            end
            
            % Calculate statistics:
            self.minV = min(values);
            self.maxV = max(values);
            self.meanV = mean(values);
            self.stdV = std(values);
        end
        
        %% Public properties:
        function Min = get.Min(self)
            Min = self.minV;
        end
        
        function Max = get.Max(self)
            Max = self.maxV;
        end
        
        function Mean = get.Mean(self)
            Mean = self.meanV;
        end
        
        function STD = get.STD(self)
            if(isnan(self.stdV))
                STD = 0;
            else
                STD = self.stdV;
            end
        end
    end
    
end

