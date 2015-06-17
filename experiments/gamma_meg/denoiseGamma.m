function denoiseGamma(subjects)
%
% denoiseGamma(subjects, howToDenoise)
%
% INPUTS:
% subjects  : Number of subjects one would like to denoise
%
% DESCRIPTION: Wrapper function to denoise multiple MEG visual gamma data
% sets. The results of the denoising are written to file.

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
epoch_start_end        = [0.050 1.049];% start and end of epoch, relative to trigger, in seconds
fs                     = 1000;        % sample rate
intertrial_trigger_num = 11;          % the MEG trigger value that corresponds to the intertrial interval
verbose                = true;

% denoise parameters (see denoisedata.m)
opt.pcchoose        = 1.05;   % denoise with exactly 10 PCs for broadband
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
opt.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis


% 2. We need a new evoked function that extracts the evoked amplitude
evokedfun             = @(x)getevoked(data,opts); % function handle to determine noise pool

% 3. We need a new eval function that extracts gamma and broadband
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband


d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

subjects = 5;

% Load and denoise data, one subject at a time
for which_subject = subjects
   
    %% Load data (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{which_subject}, 'raw'), '*Gamma*');
    
    %% Extract triggers
    trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    
    %% Make epochs
    [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
    % remove intertrial intervals
    iti               = conditions == intertrial_trigger_num;
    ts                = ts(:,~iti, :);
    conditions        = conditions(~iti);
    
    %% Remove bad epochs    
    
    [ts(:,:,data_channels), bad_channels, bad_epochs] = meg_preprocess_data(ts(:,:,data_channels), ...
        var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
    
    ts = ts(:,~bad_epochs,data_channels);
    ts(:,:, bad_epochs) = NaN;
    
    conditions = conditions(~bad_epochs);
    
    %% Make the design matrix
    num_non_blank = length(unique(conditions))-1;
    design = zeros(length(conditions), num_non_blank);
    for ii = 1:length(conditions)
        if conditions(ii) <= max(num_non_blank)
            design(ii,conditions(ii)) = 1;
        end
    end
    
    
    
    %% Denoise the data
    % Permute sensorData for denoising
    ts = permute(ts, [3 1 2]);
    
    % Denoise for broadband analysis
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,ts,evokedfun,evalfun,opt);
        
        

end
