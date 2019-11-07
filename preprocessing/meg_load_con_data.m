function [ts, meg_files] = meg_load_con_data(data_pth, fname_str)
% This function loads the CON MEG files and concatenates them, if there is
% more than one.
%
% Dependencies: fieldtrip
%
% INPUTS:
% data_pth      Absolute data path to the session folder. This folder must
%               contain a raw and behavior folder for the sqd and mat files.
% fname_str     A character array contained in the name of data files
%
% OUTPUTS:
% ts:           Timeseries of MEG experiment (timepoints x channels)
% meg_files:    MEG file names, needed to read header for topoplot

% If user doesn't specify fname_str, look for any sqd file
if ~exist('fname_str', 'var'), fname_str = '*'; end

% Find MEG and stimulus files:
meg_files = dir(fullfile(data_pth, sprintf('%s.con', fname_str)));
num_files = numel(meg_files);

% Load the MEG data for all possible blocks
ts = cell(1, num_files);
for ii = 1:num_files
    ts{ii} = ft_read_data(fullfile(data_pth, meg_files(ii).name));
end

% Concatenate timeseries if there are more blocks of MEG data
% TODO: Make this a general statement
ts = vertcat(ts{:});

return