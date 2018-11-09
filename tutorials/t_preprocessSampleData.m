%% t_preprocessSampleData

% This is a tutorial to describe general preprocessing steps starting from
% raw MEG data. Data can be downloaded from OSF website: 
%          
%       https://osf.io/rtseb/?action=download

% Overview:
%   0. Download data from OSF, move to data folder and unzip
%   1. Load data in with Fieldtrip
%   2. Get triggers
%   3. Epoch Data
%   4. (Optional) Denoise data

% Dependencies
%   * meg_utils (https://github.com/WinawerLab/meg_utils)
%   * Fieldtrip (https://github.com/fieldtrip/fieldtrip)
%   * (optional) megdenoise (https://github.com/elinekupers/denoiseproject)

% see also: t_forwardmodels

%% 0. Move downloaded data to folder and unzip

saveData   = true; % save preprocessed data later?
dataPth    = '~/Documents/';
cd(dataPth); 

% Check if data is moved, otherwise, do that and unzip compressed file 
if ~exist(fullfile(dataPth,'MEGSampleData'), 'dir')
    
    % Move zipped folder
    movefile(fullfile('~/Downloads/MEGSampleData.zip'), dataPth)

    % Go folder and unzip file
    if ~exist('MEGSampleData','dir'); mkdir('MEGSampleData'); end
    unzip(fullfile(dataPth, 'MEGSampleData.zip'), 'MEGSampleData');
    
end
    
%% 1. Load data in with Fieldtrip

dataPth = fullfile(dataPth, 'MEGSampleData');
[ts, megFiles] = meg_load_sqd_data(dataPth, '*SSMEG*'); % ts: time x channels
hdr = ft_read_header(fullfile(megFiles(1).folder, fullfile(megFiles(1).name));

%% 2. Get triggers

triggerChannels = 161:168;
trigger         = meg_fix_triggers(ts(trigger_channels,:)');


%% 3. Epoch Data

% Define length of epochs relative to trigger onset and get sample rate
epochLength = [0 1]; % start and end of epoch (seconds), zero is trigger onset
fs          = 1000;  % Hz

% Get onsets of triggers
onsets      = ssmeg_trigger_2_onsets(trigger, [], 'meg');

% Get data and conditions from raw ts
[sensorData, conditions] = meg_make_epochs(ts', onsets, epochLength, fs);

% Save preprocessed data if requested
if saveData
    if ~exist(fullfile(dataPth, 'processed'),'dir'); mkdir(fullfile(dataPth, 'processed'),'dir'); end
    save(fullfile(dataPth, 'processed', 'sensorData.mat'), 'sensorData')
    save(fullfile(dataPth, 'processed', 'conditions.mat'), 'conditions')
end

%% 4. (Optional) Denoise data

% Define frequencies of interest
f                  = 0:150;   % limit frequencies to [0 150] Hz
slFreq             = 12;      % Stimulus-locked frequency
slFreqIndex        = slFreq + 1;
tol                = 1.5;     % exclude frequencies within +/- tol of sl_freq
slDrop             = f(mod(f, slFreq) <= tol | mod(f, slFreq) > slFreq - tol);
   
% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
ln_drop            = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lf_drop            = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, abIndex]       = setdiff(f, [slDrop ln_drop lf_drop]);
keepFrequencies    = @(x) x(abIndex);
bbFrequencies      = f(abIndex);

% Get other function inputs:
design             = conditions2design(conditions);  % create design matrix (conditions x epochs) 
evokedfun          = @(x)getstimlocked(x,slFreqIndex); % function handle to determine noise pool
evalfun            = @(x)getbroadband(x,keepFrequencies,fs);  % function handle to compuite broadband with a sample rate of 1 kHz

opt.pcchoose       = -10;   % Take 10 PCs
opt.npcs2try       = [];    % if pcchoose is negative, then leave blank npcs2try, since it will only try nr of PCs specified above
opt.preprocessfun  = @(x) bbFilter(x, bbFrequencies);  % preprocess data with a filter for broadband analysis
nrControlModes     = 0;     % no controls, just the default denoising

% Denoise data now
[results, ~, ~, denoisedts] = ...
    denoisedata(design,sensorData,evokedfun,evalfun,opt);

% Save data if requested
if saveData
    if ~exist(fullfile(dataPth, 'processed'),'dir'); mkdir(fullfile(dataPth, 'processed'),'dir'); end
    save(fullfile(dataPth, 'processed', 'denoisedts.mat'), 'denoisedts')
    save(fullfile(dataPth, 'processed', 'results.mat'), 'results')
end























