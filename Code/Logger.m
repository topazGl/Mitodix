classdef Logger < handle 
    % Source: http://www.mathworks.com/matlabcentral/fileexchange/33532-log4matlab
    % Copyright (c) 2013, Gavin
    % All rights reserved.     
    properties(SetObservable = true)
        mlog;
    end
    
    properties(Access = protected)
        logger;          
    end
    
    properties(Access = private)
        listener2mlog;
    end
    
    properties(SetAccess = protected)
        fullpath = 'Mitodix.log';
        append = true;
        % By default nothing printed to file will be printed to the console. Must be greater or equal to log to file level.
        commandWindowLevel = LogType.NONE; 
    end

    
    methods
%% Constructor
        function self = Logger(fullpath_passed,doAppend)           
            if 0 < nargin
                self.SetFilename(fullpath_passed);
                if nargin > 1
                    self.SetAppend(doAppend);
                end
            end
            self.listener2mlog = addlistener(self,'mlog','PostSet',@self.WriteToFile);
        end
 
%% Public properties:
        function SetFilename(self,fullpath_passed)
        % set fileName
            self.fullpath = fullpath_passed;
            % If not appending, create a new log file:
            if ~self.append
                fid = fopen(self.fullpath,'w');
                fclose(fid);
            end
        end
  
        function SetLoggerLevel(self,loggerIdentifier,value)%#ok<MANU>
        % set loggerLevel:
            loggerIdentifier = strrep(loggerIdentifier, '.', '_');
            eval(['self.logger.',loggerIdentifier,'.value = ',num2str(value),';']);
            pause(0);
        end
              
        function result = GetLoggerLevel(self,loggerIdentifier)
        % get loggerLevel:
            % Default is debug:
            result = LogType.DEBUG;            
            loggerIdentifier = strrep(loggerIdentifier, '.', '_');
            
            if isfield(self.logger,loggerIdentifier) 
                eval(['result = self.logger.',loggerIdentifier,'.value;']);
            else
                self.SetLoggerLevel(loggerIdentifier,result)
            end                
        end
        
        function SetCommandWindowLevel(self,value)
        % set commandWindowLevel:
            self.commandWindowLevel = value;
            pause(0);
        end
    end
 %% Private methods:    
    methods (Access = private)

        function SetAppend(self,value)
        % set append
            if value
                self.append = true;
            else
                self.append = false;
            end
        end
               
        function WriteToFile(self,src,evnt)             %#ok<INUSD>
        % write to file:    
            if ~iscell(self.mlog) || size(self.mlog,2)~=3
                error('Problem with the value of mlog that has been set');
            elseif isempty(self.fullpath)
                error('You must set fullpath using SetFilename() method');
            end
            
            loggerIdentifier = strrep(self.mlog{2}, '.', '_');            
            % If there is no field by this name or the logger level set is 
            % less than the level passed then print to file
            DoIt = true;
            if isfield(self.logger,loggerIdentifier) 
                tempField = getfield(self.logger,loggerIdentifier); %#ok<GFLD>
                if self.mlog{1} < tempField.value
                    DoIt = false;
                end
            end
                
            if DoIt
                if (self.mlog{1} == LogType.DEBUG), levelStr = 'DEBUG';
                elseif (self.mlog{1} == LogType.WARNING), levelStr = 'WARNING';
                else levelStr = 'ERROR';
                end

                try                 
                    fid = fopen(self.fullpath,'a');
                    fprintf(fid,'%s %s %s - %s\r\n' ...
                        , datestr(now,'yyyy-mm-dd HH:MM:SS.FFF') ...
                        , levelStr ...
                        , self.mlog{2} ... % Have left this one with the '.' if it is passed
                        , self.mlog{3});
                    fclose(fid);
                catch ME_1
                    display(ME_1);
                end  
				
                % If necessary write to command window (note it will only be written to the conole if: 
                % 1) it is written to file AND
                % 2) the min commandWindowLevel is less than or equal to the log level
                if self.commandWindowLevel <= self.mlog{1}
                    display([self.mlog{2},': ',self.mlog{3}]);                
                end
            end

        end
    end

    %% Static methods:
	methods (Static)
				
		function str = MatrixToString(matrix)
        % matrix to string:
			str = char(10); % Line feed
			for (i=1:size(matrix,1))
				for (j=1:size(matrix,2))
					element = matrix(i,j);
					str = sprintf([str,'% 10.4f   '],element);
				end
				str = [str,char(10)]; %#ok<AGROW>
			end		
		end
		
		function str = ExceptionToString(ME)
        % exception to string:
			str = char(10); % Line feed			
			str = [str,'ME.identifier: ',ME.identifier,char(10)];
			str = [str,'ME.message: ',ME.message,char(10)];
			for i=1:length(ME.stack)
				str = [str,'ME.stack: ',ME.stack(i).name,', Line ',num2str(ME.stack(i).line),char(10)]; %#ok<AGROW>
			end
        end
        
        function str = Me()
        %> Returns the current function name
            dbstackOutput = dbstack;
            if size(dbstackOutput,1) < 2
                str = '';
                return 
            end
            fullName = dbstackOutput(end).name;
            strReversed = strtok(fullName(end:-1:1),'.');
            str = [strReversed(end:-1:1),': '];
        end
        
    end	
end


