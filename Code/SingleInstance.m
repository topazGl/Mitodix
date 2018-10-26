classdef (Sealed) SingleInstance < handle
    % SINGLEINSTANCE this is based directly on the singleinstance class
    % from matlab documentation. Look up singleinstance for more details
    % Copyright (c) 2013, Gavin
    % All rights reserved.
    properties (Access = protected)
        % This is the general value that will be stored for whatever the 
        value;
        requiredType;
    end
    methods
        function delete(self)
            delete(self.value);
        end
    end
    methods (Access = private)
        function self = SingleInstance(requiredType_passed)
            self.requiredType = requiredType_passed;
        end
        
%% testAndSet: Handles the general input
        function testAndSet(self,varagin)
            if isa(varagin,self.requiredType)
                if(isequal('Log', self.requiredType))
                    self.value = varagin;
                end
            else
                whosOutput = whos('varagin');
                warning(['Expecting varagin to be of type: ',self.requiredType ...
                        ,'. However, input passed in is of type: ',whosOutput.class]); %#ok<WNTAG>
            end
            
            if(isa(self.value, 'Config'))
                if (isempty(self.value) || ~isequal(self.value.DataName, varagin))
                    self.value = Config(varagin);
                end
            end
        end        
    end
    
    methods (Static)
		
%% Logger 
        function result = Logger(varagin)
            persistent localObj;
            if (isempty(localObj) || ~isvalid(localObj) || isempty(localObj.value)  || ~isvalid(localObj.value))
                localObj = SingleInstance('Log'); 
                localObj.value = Logger('Mitodix.log');
            end
            
            if nargin > 0
                localObj.testAndSet(varagin);
            end
            result = localObj.value;
        end	
        
%% Config 
        function result = Config(varagin)
            persistent localConfig;
            if (isempty(localConfig) || ~isvalid(localConfig) || isempty(localConfig.value) || ~isvalid(localConfig.value))
                localConfig = SingleInstance('char');
                localConfig.value = Config.empty;   
            end
            
            if nargin > 0
                localConfig.testAndSet(varagin);
            end
            result = localConfig.value;
		end	
    end
end