%% t_preprocessSampleData


% This is a tutorial to describe general preprocessing steps starting from
% raw MEG data.

% Overview:
%   0. Download data from OSF
%   1. Load data in with Fieldtrip
%   2. Get triggers
%   3. Epoch Data
%   4. (Optional) Denoise data
%
%   To add dataset to Brainstorm and get gain matrix: 
%   x. 

% Dependencies
%   * meg_utils (https://github.com/WinawerLab/meg_utils)
%   * Fieldtrip (https://github.com/fieldtrip/fieldtrip)

%   * (optional) megdenoise (https://github.com/elinekupers/denoiseproject)

%% 0. Download data from OSF

% Site to retrieve the data
dirProject = 'https://osf.io';
urlStr     = 'rtseb';

% Folder to save data
writePth    = '~/Downloads/';
readPth    = fullfile(dirProject, urlStr, '?action=download');

% Download data
websave(writePth,readPth);

% Move and unzip folder
dataPth    = '~/Documents/';
mv(fullfile(writePth, 'MEGSampleData.zip'), dataPth)
unzip(fullfile(dataPth, 'MEGSampleData.zip'));

%% 1. Load data in with Fieldtrip

dataPth = fullfile(writePth, 'SSMEG');
[ts, megFiles] = meg_load_sqd_data(dataPth, '*SSMEG*'); % ts: time x channels


%% 2. Get triggers

triggerChannels = 161:168;
trigger = meg_fix_triggers(ts(trigger_channels,:)');

%% 3. Epoch Data

epochLength = [0 1]; % start and end of epoch (seconds), zero is trigger onset
fs          = 1000;  % Hz

% Get onsets of triggers
onsets      = ssmeg_trigger_2_onsets(trigger, [], 'meg');

% Get data and conditions from raw ts
[sensorData, conditions] = meg_make_epochs(ts', onsets, epochLength, fs);

%% 4. (Optional) Denoise data

























