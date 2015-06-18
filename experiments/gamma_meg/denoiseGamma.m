function denoiseGamma(subject)
%
% denoiseGamma(subject)
%
% INPUTS:
% subjects  : Subject number to denoise
%
% DESCRIPTION: Wrapper function to denoise  MEG visual gamma data sets.


% TODO:
%      1. function should return something
%      2. write a preprocessing function
%      3. add environmental denoising

project_pth                     = '/Volumes/server/Projects/MEG/Gamma/Data';

% data to be analysed
data_pth                      = '*_Gamma_*subj*';


% preprocessing parameters (see dfdPreprocessData)
var_threshold          = [0.05 20];
bad_channel_threshold  = 0.2;
bad_epoch_threshold    = 0.2;
data_channels          = 1:157;
environmental_channels = 158:160;
trigger_channels       = 161:164;
epoch_start_end        = [0.050 1.049-.5];% start and end of epoch, relative to trigger, in seconds
fs                     = 1000;        % sample rate
intertrial_trigger_num = 11;          % the MEG trigger value that corresponds to the intertrial interval
blank_condition        = 10;          % the MEG trigger value that corresponds to trials with zero contrast
verbose                = true;



d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

% Load and denoise data, one subject at a time

%% Load data (SLOW)
raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject}, 'raw'), '*Gamma*');

%% Extract triggers
trigger = meg_fix_triggers(raw_ts(:,trigger_channels));

%% Make epochs
[ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
% remove intertrial intervals
iti               = conditions == intertrial_trigger_num;
ts                = ts(:,~iti, :);
conditions        = conditions(~iti);
conditions(blank_condition) = 0;

%% Remove bad epochs and bad channels

[ts(:,:,data_channels), bad_channels, bad_epochs] = meg_preprocess_data(ts(:,:,data_channels), ...
    var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);

ts = ts(:,~bad_epochs,~bad_channels);

conditions = conditions(~bad_epochs);

design_mtrx = conditions2design(conditions);

%% Denoise the data

% denoise parameters (see denoisedata.m)
opt.pcchoose        = 1.05;
opt.npcs2try        = 10;
opt.resampling      = {'xval','xval'}; % could be {'boot' 'boot'};
opt.pcselmethod     = 'r2';            % could be 'snr';

%% TODO

% 1. Change the preprocess function
%        (a) We want to remove the first few hundred ms (maybe 300?)
%        (b) We want to keep these frequencies and filter out the rest:
%               f_use4fit = f((f>=35 & f < 40) |(f > 40 & f <= 57) | ...
%                   (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
%
t                   = epoch_start_end(1):1/fs:epoch_start_end(2);
f                   = (0:length(t)-1)*fs/length(t);
keep_timepts        = t(t>.250);
keep_frequencies    = f((f>=35 & f < 40) |(f > 40 & f <= 57) | ...
                   (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
opt.preprocessfun   = @(x)gammapreprocess(x, keep_timepts, keep_frequencies);  % preprocess data by clipping 1st

% 2. We need a new evoked function that extracts the evoked amplitude
evokedfun           = @(x)getevoked(x, fs, design_mtrx, [-30 30]); % function handle to determine noise pool

% 3. We need a new eval function that extracts gamma and broadband
f                     = (0:size(ts,1)-1)/max(t);
freq                  = megGetSLandABfrequencies(f,[],60);
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband


% Permute sensorData for denoising
ts = permute(ts, [3 1 2]);

% Denoise for broadband analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design_mtrx,ts,evokedfun,evalfun,opt);




