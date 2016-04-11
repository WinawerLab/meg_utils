function [ts, conditions, badChannels, badEpochs] = gamma_preprocess_raw(sessionNum)
% [ts, conditionsVector] = gamma_preprosess_raw(sessionNum)
% 
% {For use in Nicholas' new pipeline}
% loads the squid file
% fixes trigger anomalies
% epoches time series into a 3D (time x epoch x channel) matrix
% 
% Input : sessionNum of the data being analysed
% Output: ts - 3D time series (time x epoch x channel)
%       : conditions - vector of size(ts,2) with trigger values
%       : badChannels - binary vector corresponding to channels that were
%       removed
%       : badEpochs - binary vector representing epochs that were removed

%% preprocessing variables

data_channels     = 1:157;
trigger_channels  = 161:164;
fs                = 1000; % sampling rate in milliseconds
epoch_start_end   = [0.050 1.049]; % start and end of epoch in seconds

% parameters for the removal of nuissance data
var_threshold         = [0.05 20]; 
bad_channel_threshold = 0.2;
bad_epoch_threshold   = 0.2;
verbose               = true;

SAVE_PREPROCESSED = false;

%% Get session paths

% parent dir to this session's folder on the server
sessionPath = meg_gamma_get_path(sessionNum);

%% Load sqd file

raw = meg_load_sqd_data(fullfile(sessionPath, 'raw'), '*Gamma*');

%% Extract trigger sequence

trigger = meg_fix_triggers(raw(:,trigger_channels));

%% Partition into epochs

[ts, conditions] = meg_make_epochs(raw, trigger, epoch_start_end, fs);

%% Find and return bad channels/epochs

[ts(:,:,data_channels), badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
    var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);


%% Save if needed 

if SAVE_PREPROCESSED
    savePath    = fullfile(sessionPath, 'processed');
    if ~exist(savePath, 'dir'), mkdir(savePath); end
    fileName = fullfile(savePath, sprintf('s%03d_preprocessed_ts_%s.mat', sessionNum, datestr(now, 'mm.dd.yy')));
    save(fileName, 'ts', 'conditions', 'badChannels', 'badEpochs');
end





end