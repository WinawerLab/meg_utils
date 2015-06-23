%% denoiseGamma scripts
%
% DESCRIPTION: Wrapper script to denoise  MEG visual gamma data sets.


% ------------------------------------------------------------------------
% ----------------- Check options & define variables----------------------
% ------------------------------------------------------------------------

% Where to find data?
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';

% Type of data
data_pth                      = '*_Gamma_*subj*';

% Subject number to analyze
subject                       = 5; % Fifth folder in project path, not neccarily session 5.

% preprocessing parameters (see dfdPreprocessData)
var_threshold                 = [0.05 20];
bad_channel_threshold         = 0.2;
bad_epoch_threshold           = 0.2;
data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;
epoch_start_end               = [0.050 1.049];% start and end of epoch, relative to trigger, in seconds
fs                            = 1000;        % sample rate
intertrial_trigger_num        = 11;          % the MEG trigger value that corresponds to the intertrial interval
blank_condition               = 10;          % the MEG trigger value that corresponds to trials with zero contrast
verbose                       = true;
denoise_with_nonphys_channels = false;

% Find subject path
d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

% ------------------------------------------------------------------------
% --------------------- Load & preprocess data ---------------------------
% ------------------------------------------------------------------------


%% Load data (SLOW)
raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject}, 'raw'), '*Gamma*');

%% Extract triggers
trigger = meg_fix_triggers(raw_ts(:,trigger_channels));

%% Make epochs
%   ts is [time points by epochs x channels]
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


% Denoise data with 3 noise channels
if denoise_with_nonphys_channels
    ts = meg_environmental_denoising(ts, environmental_channels,data_channels(~bad_channels));
end


%% -----------------------------------------------------------------------
% -------------------------- Denoise data --------------------------------
% ------------------------------------------------------------------------

% Define denoise Parameters (see denoisedata.m)
opt.pcchoose        = 1.05;
opt.npcs2try        = 10;
opt.resampling      = {'xval','xval'}; % could be {'boot' 'boot'};
opt.pcselmethod     = 'r2';            % could be 'snr';
if verbose; opt.verbose = true; end

t                   = epoch_start_end(1):1/fs:epoch_start_end(2);
f                   = (0:length(t)-1)*fs/length(t);
keep_timepts        = t(t>.250);
keep_frequencies    = f((f>=35 & f < 40) |(f > 40 & f <= 57) | ...
                   (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
opt.preprocessfun   = @(x)gammapreprocess(x, t, f, keep_timepts, keep_frequencies);  % preprocess data by clipping 1st

% A new evoked function that extracts the evoked amplitude
evokedfun           = @(x)getevoked(x, fs, design_mtrx, [-30 30]); % function handle to determine noise pool

% Eval function that extracts broadband
f                   = (0:length(keep_timepts)-3)*fs/length(keep_timepts);
freq                = megGetSLandABfrequencies(f,[],60);
evalfun             = @(x)getbroadband(x,freq);  % function handle to compuite broadband

% Permute sensorData for denoising
%   [time points by epochs x channels] => [channels x time points x epochs]
ts = permute(ts, [3 1 2]);

% Denoise for broadband analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design_mtrx,ts,evokedfun,evalfun,opt);



