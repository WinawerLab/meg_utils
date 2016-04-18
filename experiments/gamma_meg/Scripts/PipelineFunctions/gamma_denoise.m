function [tsDenoised] = gamma_denoise(ts, conditions, sessionNum)
% [denoisedTs] = gamma_denoise(ts, conditions, sessionNum, badChannels, badEpochs)
% wrapper function for Eline's PCA denoising program
% takes the epoched time series and returns a denoised version of the same
% format
% additional information about the denoising process can be saved to file
%
% Inputs  : ts - epoched time series that has been preprocessed
%         : condVector - a vector of trigger values corresponding to the
%        epochs in the ts. condVector must be size(ts, 2)
% Outputs : denoisedTs - epoched time series that has been denoised

%% Options and parameters

EPOCH_START_END = [0.050 1.049];
fs              = 1000;
subject         = sessionNum;

if subject >= 17
    ITI       = 14;
else
    ITI       = 10;
end
blankCond = ITI -1;

save_data = true;
verbose = true;

%% get paths
% path to denoising project
addpath(genpath('~/matlab/git/denoiseproject'))
% path so session folder
sessionPath = meg_gamma_get_path(sessionNum);
% path to filed trip
if isempty(which('sqdread')),
    meg_add_fieldtrip_paths(ft_pth,{'yokogawa', 'sqdproject'})
end

%% Get time series ready for denoising
% remove intertrial epochs
ITIepochs  = conditions == ITI;
ts         = ts(:, ~ITIepochs, :);
conditions = conditions(conditions ~= ITI);

% change the trigger value of blank stim to 0
conditions(conditions == blankCond) = 0;
%% Make a design matrix

design_mtrx = conditions2design(conditions);


%% -----------------------------------------------------------------------
% -------------------------- Denoise data --------------------------------
% ------------------------------------------------------------------------

% Define denoise Parameters (see denoisedata.m)
opt.pcchoose        = 1.05;
opt.npcs2try        = 10;
opt.resampling      = {'xval','xval'}; % could be {'boot' 'boot'};
opt.pcselmethod     = 'r2';            % could be 'snr';
if verbose; opt.verbose = true; end

% time and frequencies for full-length epochs (prior to truncating)
evoked_cutoff       = .250; % use only time points >= evoked_cutoff
t                   = EPOCH_START_END(1):1/fs:EPOCH_START_END(2);
t_clipped           = t(t >= evoked_cutoff);
f                   = (0:length(t_clipped)-1)*fs/length(t_clipped);
%     keep_frequencies    = @(x) x((x>=35 & x < 40) |(x > 40 & x <= 57) | ...
%         (x>=65 & x <= 115) | (x>=126 & x<= 175) | (x>=186 & x<= 200));

keep_frequencies = @(x) x((x>=35 & x <= 57) | (x>=63 & x <= 115) | (x>=126 & x <= 175) | (x>=186 & x <= 200));

% preprocess data by clipping 1st n time points (to eliminate evoked
% response) and by removing all frequencies that will not be used to fit
% gamma and broadband later on.
opt.preprocessfun   = @(x)gammapreprocess(x, t, f, evoked_cutoff, keep_frequencies(f));

% Evoked function is used to identify the noise pool. It finds the peak
% evoked signal in each epoch for each channel.
evokedfun           = @(x)getevoked(x, fs, design_mtrx, [-30 30]); % function handle to determine noise pool

% Eval function that extracts broadband
evalfun             = @(x)getbroadband(x,keep_frequencies, fs);  % function handle to compute broadband

% Permute sensorData for denoising
%   [time points by epochs x channels] => [channels x time points x epochs]
ts = permute(ts, [3 1 2]);

% Denoise for broadband analysis
[results,evalout,denoisedspec,tsDenoised] = denoisedata(design_mtrx,ts,evokedfun,evalfun,opt);

if save_data
    thisDate = datestr(now, 'mm.dd.yy');
    fileName = sprintf('s_%03d_denoisedData%s.mat', subject, thisDate);
    save(fullfile(meg_gamma_get_path(subject), 'processed', fileName),...
        'results','evalout','tsDenoised','opt', denoisedspec)
end
end