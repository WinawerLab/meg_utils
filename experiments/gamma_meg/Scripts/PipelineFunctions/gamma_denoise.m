function [tsDenoised] = gamma_denoise(ts, opt)
% 
% [tsDenoised] = gamma_denoise(ts, opt)
% Function that calls MEG denoise with the standard parameters and denoises
% up to 10 pcs. Function outputs denoised timeseries only, no results from
% evaluation function, because this contains both broadband and gamma. A
% separate modelfit function (gamma_spectral_analysis.m) is necessary.
%
% INPUTS:
%       ts    :    3D array with epoched time series that has been
%                  preprocessed (time points x epochs x channels)
%       opt   :    structure that contains used experimental/analysis
%                   parameters
%                                
% OUTPUTS:
%       denoisedTs : 3D array with epoched time series that has been
%                  denoised (time points x epochs x channels)

% First version by Nicholas Chua (April 2016)
%       07.12.16 Clean up (EK)

%% Get time series ready for denoising

% Change the trigger value of blank stim to 0
opt.params.conditions(opt.params.conditions == opt.params.baselineCondition) = 0;

%% Make a design matrix
designMtrx = conditions2design(opt.params.conditions);

%% -----------------------------------------------------------------------
% -------------------------- Denoise data --------------------------------
% ------------------------------------------------------------------------

% Define denoise Parameters
opt.pcchoose        = -10;
opt.npcs2try        = [];
opt.resampling      = {'xval','xval'}; % could be {'boot' 'boot'};
opt.pcselmethod     = 'r2';            % could be 'snr';

% Truncate time points in epochs since we use first 250 ms to get evoked
% field in order to define the channels that go into the noise pool
opt.evokedCutOff    = .250; % use only time points >= evoked_cutoff
t                   = opt.params.epochStartEnd(1):(1/opt.fs):opt.params.epochStartEnd(2);
tClipped            = t(t >= opt.evokedCutOff);
f                   = (0:length(tClipped)-1)*opt.fs/length(tClipped);
keepFrequencies     = @(x) x((x>=35 & x <= 56) | (x>=63 & x <= 115) | (x>=126 & x <= 175) | (x>=186 & x <= 200)); % Todo: extract frequencies from opt.fitFreq so we don't have to define it again here

% Preprocessing function will clip time points up to opt.evokedCutOff
% (to eliminate evoked response) and by removing all frequencies that will 
% not be used to fit gamma and broadband later on.
opt.preprocessfun   = @(x)gammapreprocess(x, t, f, opt.evokedCutOff, keepFrequencies(f));

% Evoked function will identify the noise pool, by finding the peak evoked 
% signal in each epoch for each channel.
evokedfun           = @(x)getevoked(x, opt.fs, designMtrx, [-30 30]);

% Evaluation function that extracts broadband (and most likely some narrow-band gamma)
evalfun             = @(x)getbroadband(x, keepFrequencies, opt.fs);

% Permute sensorData for denoising
ts = permute(ts, [3 1 2]);

% Denoise for broadband analysis
[results,~,~,tsDenoised] = denoisedata(designMtrx,ts,evokedfun,evalfun,opt);

% Set things right before we move on
if opt.MEGDenoise==0; opt.MEGDenoise = 1; end
tsDenoised = permute(tsDenoised{1}, [2 3 1]); % back to time x epochs x chan

if opt.saveData
    thisDate = datestr(now, 'mm.dd.yy');
    fileName = sprintf('s_%02d_tsDenoised_%s.mat', sessionNum, thisDate);
    save(fullfile(meg_gamma_get_path(sessionNum), 'processed', fileName),...
        'results','tsDenoised','opt')
end
end