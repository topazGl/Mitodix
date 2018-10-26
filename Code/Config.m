classdef Config < handle
% Configuration class

    properties (SetAccess = protected)
        DataName;               % data set name / identification.    
        InputImPath;            % input gray-level directory path.
        InputImExtension;       % input gray-level files extension.
        InputImExpression;      % input gray-level file name expression.
        InputSegPath;           % input segmentation directory path.
        InputSegExtension;      % input segmentation files extension.
        InputSegExpression;     % input segmentation file name expression.
        OutputPath;             % output directory name
        GTPath;                 % Ground truth *.mat path
        MinCellLength;          % min cell major axis (in pixels).
        FrameRate;              % seauence frame rates [frame per minutes].
        CellCycle;              % Cell cycle length [frames].
    end
    
    properties (Access = protected)        
        cellCycleTime;        % Cell cycle length [hours].          
    end

    methods      
        % Constructor
        function self = Config(dataName, downscaleFactor, dt)
        % Constructor
            if(~exist('dataName', 'var') || isempty(dataName))
                xmlFile = xmlread('Mitodix.config');
                factory = javax.xml.xpath.XPathFactory.newInstance();
                xpath = factory.newXPath();
                expr = xpath.compile('root//chosen/@runName');
                val = expr.evaluate(xmlFile);
                dataName = char(val);
            end    
            if(~exist('downscaleFactor', 'var') || isempty(downscaleFactor))
                downscaleFactor = 1;
            end
            
            if(~exist('dt', 'var') || isempty(dt))
                dt = 1;
            end
            
            self.DataName = dataName;
            self.ChooseSection(self.DataName, downscaleFactor, dt);
        end

        % Perform subclass save operations
        function sobj = saveobj(self)
            sobj = struct(self);          
        end
  
        % Get object as struct
        function sobj = struct(self)
            sobj = [];
            sobj.DataName = self.DataName;
            sobj.InputImPath = self.InputImPath;            
            sobj.InputImExtension = self.InputImExtension;       
            sobj.InputSegPath = self.InputSegPath;           
            sobj.InputSegExtension = self.InputSegExtension;      
            sobj.OutputPath = self.OutputPath;             
            sobj.cellCycleTime = self.cellCycleTime;     
            sobj.FrameRate = self.FrameRate;         
            sobj.CellCycle = self.CellCycle;       
            sobj.GTPath = self.GTPath;      
            sobj.MinCellLength = self.MinCellLength;
        end
        
        % Choose config xml section according to chosen data set ("run" tag)
        function ChooseSection(self, dataName, downscaleFactor, dt)
        % Choose config xml section according to chosen data set ("run" tag)
        %   ChooseSection(self, dataName)
        % Inputs:
        %   self = Config object.
        %   dataName = data set section name in config xml file ("run" tag).
            self.DataName = dataName;
            try
                % Read configuration XML file:
                xmlFile = xmlread('Mitodix.config');
                factory = javax.xml.xpath.XPathFactory.newInstance();
                xpath = factory.newXPath();
                
                % Get XML attributes:
                self.InputImPath = fullfile(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputImage', 'path')));
                self.InputSegPath = fullfile(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputSegmentation', 'path')));
                self.OutputPath = fullfile(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'output', 'path')));
                self.GTPath = fullfile(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'groundTruth', 'path')));

                self.cellCycleTime = str2double(Config.GetAttributeValue(xmlFile, xpath, dataName, 'cellCycleTime', 'hours'));
                if (~isempty(self.cellCycleTime))
                    self.cellCycleTime = 24;
                end
                self.FrameRate = str2double(Config.GetAttributeValue(xmlFile, xpath, dataName, 'frameRate', 'minutes'));
                if (isempty(self.FrameRate))
                    self.FrameRate = 20;
                end
                self.FrameRate = self.FrameRate * dt;
                self.CellCycle = (60*self.cellCycleTime)/self.FrameRate;
                
                self.MinCellLength = str2double(Config.GetAttributeValue(xmlFile, xpath, dataName, 'minCellLength', 'pixels'));
                self.MinCellLength = floor(self.MinCellLength / downscaleFactor);
                
                inputExt = lower(strtrim(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputImage', 'extension'))));
                if ((~isempty(inputExt)) && (inputExt(1) == '.'))
                    self.InputImExtension = inputExt(2:end);
                else
                    self.InputImExtension = inputExt;
                end
                
                self.InputImExpression = char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputImage', 'expression'));
                
                inputExt = lower(strtrim(char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputSegmentation', 'extension'))));
                if ((~isempty(inputExt)) && (inputExt(1) == '.'))
                    self.InputSegExtension = inputExt(2:end);
                else
                    self.InputSegExtension = inputExt;
                end
                
                self.InputSegExpression = char(Config.GetAttributeValue(xmlFile, xpath, dataName, 'inputSegmentation', 'expression'));
            catch EX
                throw(EX);
            end
        end        
    end
    
    methods(Static)
        % Load config object from file
        function obj = loadobj(objStruct)
            obj = Config(objStruct.DataName);
        end
        
        % Get xml attribute value
        function val = GetAttributeValue(xmlFile, xpath, dataName, tag, field)
        % Get xml attribute value
        %   val = GetAttributeValue(xmlFile, xpath, dataName, tag, field)
        % Inputs:   xmlFile = xml file object (DeferredDocumentImpl).
        %           xpath = XPathEvaluator.
        %           dataName = the data set name ("run" tag).
        %           tag = XML tag name.
        %           field = XML tag attribute name.
        % Output:   XML attribute value for the dataName of the
        %           specific Config object.
        	attribute = [tag, '/@', field];
        	query = sprintf('root//run[@name=''%s'']/%s', dataName, attribute);
        	expr = xpath.compile(query);
        	val = expr.evaluate(xmlFile);
        end
    end

end

