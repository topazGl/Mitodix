%%%
%	This is the main script for the Mitosis detection presented in:
%	Gilad et al, Fully Unsupervised Symmetry-Based Mitosis Detection in Time-Lapse Cell Microscopy, 2018.
% (c) T.Gilad, 2018.
%
% Make sure the Backup_Mitodix.config file is updated and contains
% a section of your data name with all relevane paths and file formats.
% Ground truth path is only required in debug mode.
% cellCycleTime and frameRate - are only used to determine the frame rate
% of the output video (*.avi).
%
% The algorithm detect symmetric cell divisions, provided that the input
% segmentation is successfull in identifing single cells and separates clustered cells.
%%%
close all;
addpath(genpath('.'));
clc;

DEBUG = 0;
startFrame = []; endFrame = []; % for empty values all range is read from folder.

dataSets = {'SIM+_2D_01', 'SIM+_2D_02', 'SIM+_3D_01', 'SIM+_3D_02'};
L = length(dataSets);
comperToSvd = 0;
downscaleFactor = 1;	% change to dt > 1 if you wish to downsample data (improved runtime). 
dt = 1;					% cahnge to dt > 1 (integer) if you wish to sample the time sequance (improved runtime).
for ll=1:L
    dataName = dataSets{ll};
    fprintf(['\n***********  ', dataName, ', downscale ', num2str(downscaleFactor), ', dt=', num2str(dt), ':  ***********\n']);
    try
        copyfile('Backup_Mitodix.config', 'Mitodix.config');
        mtdx = Mitodix(DEBUG, dataName, downscaleFactor);
        mtdx.Detection(comperToSvd, startFrame, endFrame);
    catch EX
        disp(EX.message);
    end
    close all;
end