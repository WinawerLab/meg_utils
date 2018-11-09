%% t_preprocessSampleData


% This is a tutorial to describe general preprocessing steps starting from
% raw MEG data. Data can be downloaded from OSF website: 
%          
%       https://osf.io/rtseb/?action=download

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

%% 0. Move downloaded data to folder and unzip

% Move and unzip folder
dataPth    = '~/Documents/';
movefile(fullfile('~/Downloads/MEGSampleData.zip'), dataPth)

% Go folder and unzip file
cd(dataPth); if ~exist('MEGSampleData','dir'); mkdir('MEGSampleData'); end
unzip(fullfile(dataPth, 'MEGSampleData.zip'), 'MEGSampleData');

%% 1. Load data in with Fieldtrip

dataPth = fullfile(dataPth, 'MEGSampleData');
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

























