function [ts, meg_files, stim_file, num_runs] = ssmeg_loaddata(data_pth)

% This function loads the SQD MEG and Matlab stimulus files.

% Make sure that the folderstructure has a 'raw' and a 'behavior' folder. 
% Raw contains the sqd files with ssmeg in the filename. 
% Behavior contains the matlab stimulus files.

% You need the fieldtrip toolbox to run this function

% INPUTS:
% data_pth      Absolute data path to the session folder. This folder must
%               contain a raw and behavior folder for the sqd and mat files.

% OUTPUTS:
% ts:           Timeseries of MEG experiment (timepoints x channels)
% meg_files:    MEG file names, needed to read header for topoplot
% stim_file:    Matlab stimulus file of the first run. Contains parameters of
%               the experiment.
% num_runs:     Number of runs in the whole experiment.

% Go to data path
cd(data_pth)

% Find MEG and stimulus files:
if exist('/raw/*_SSMEG_Block*.sqd','dir');
        meg_files = dir('raw/*_SSMEG_Block*.sqd');
else meg_files = dir('raw/*_SSMEG_*.sqd');
end

stim_files = dir('behavior/*.mat');
num_runs   = numel(stim_files);

% Load the MEG data for all possible blocks
for nr_files = 1:numel(meg_files)
    ts{nr_files} = sqdread(fullfile('raw', meg_files(nr_files).name));
end

% Concatenate timeseries if there are more blocks of MEG data
% TODO: Make this a general statement
if numel(ts) > 1
    ii = 1:numel(ts);
%     ts = vertcat(ts{ii(1)},ts{ii(2)});
    ts = vertcat(ts{ii(:)});
else
    ii = 1:numel(ts);
    ts = vertcat(ts{ii(1)});
end

% Load one stim file, not sure whether we need all stim_files
stim_file = load(fullfile('behavior',stim_files(1).name));
