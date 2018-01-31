function [ts, opt] = gamma_preprocess_raw(sessionNum, opt)
% [ts, opt] = gamma_preprocess_raw(sessionNum, opt)
% 
% {For use in Nicholas' new pipeline}
% loads the squid file
% fixes trigger anomalies
% epoches time series into a 3D (time x epoch x channel) array
% 
% INPUTS 
%    sessionNum  :    number of dataset being analysed
%    opt         :    options for preprocessing

% OUTPUTS
%    ts          :    3D time series (time x epochs x channels)
%    opt struct added with:
%       conditions  :    vector of size(ts,2) with trigger values
%       badChannels :    binary vector representing channels that were removed
%       badEpochs   :    binary vector representing epochs that were removed

%% Load sqd file
ts = meg_load_sqd_data(fullfile(opt.sessionPath, 'raw'), '*Gamma*');

% Check dimensions ts (assuming that there are more timepoints than
% channels in a dataset
if size(ts,1) < size(ts,2)
    ts = ts';
end
%% Extract trigger sequence
trigger = meg_fix_triggers(ts(:,opt.triggerChannels));

%% Partition into epochs
[ts, conditions] = meg_make_epochs(ts, trigger, opt.params.epochStartEnd, opt.fs);

% Denoise with three environmental channels
if opt.environmentalDenoising
    ts = meg_environmental_denoising(ts);
end

%% Find and return bad channels/epochs
[ts, badChannels, badEpochs] = dfdPreprocessData(ts(:,:,opt.dataChannels), ...
    opt.varThreshold, opt.badChannelThreshold, opt.badEpochThreshold, opt.verbose);

opt.params.conditions  = conditions;
opt.params.badChannels = badChannels;
opt.params.badEpochs   = badEpochs;

%% Save if needed 
if opt.saveData
    savePath    = fullfile(opt.sessionPath, 'processed');
    if ~exist(savePath, 'dir'), mkdir(savePath); end
    fileName = fullfile(savePath, sprintf('s%03d_preprocessed_ts_%s.mat', sessionNum, datestr(now, 'mm.dd.yy')));
    save(fileName, 'ts', 'opt', '-v7.3');
end

end