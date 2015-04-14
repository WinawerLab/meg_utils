%% Script to analyze EEG pilot data

% Description comes here: Script to analyze a pilot steady state EEG data
% set (subject wl_s004). Contrast patterns were contrast-reversed at 12 Hz
% in 6 second blocks alternating with 6 s of blank (mean luminance) while
% subjects fixated the middle of the screen and detected a fixation color
% change. Stimuli consisted of either full-field (11 deg radius?? check
% this), left field, or right field apertures.
%
%
% Some definitions:
%   a run:    a period of 72 seconds during which we run one experiment from
%               vistadisp
%   a block:  a 6 second period in which the stimulus condition is the same
%   an epoch: a one second periond in which the data are binned for
%               analysis
% There are 12 blocks per run, and 6 epochs per block. There are 12 images
% per epoch because the stimulus contrast reverses 12 times per second (= 6
% Hz f1)
%
% Dependencies: meg_utils github repository


%% Define variables for this experiment and analysis
project_path          = '/Volumes/server/Projects/EEG/SSEEG/';
s_rate_eeg            = 1000;     % sample rate of the eeg in Hz
s_rate_monitor        = 60;       % sample rate of the monitor in Hz
plot_figures          = true;     % Plot debug figures or not?
images_per_block      = 72;       % number of images shown within each 6-s block of experiment
epochs_per_block      = 6;        % bin data into 1 second epochs (blocks are 6 s)
blocks_per_run        = 12;       % number of blocks in one experimental run
var_threshold         = [.05 20]; % acceptable limits for variance in an epoch, relative to median of all epochs
bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
data_channels         = 1:128;
verbose               = true;


%% Define variables for this particular subject's session
session_name   = 'SSEEG_20150403_wl_subj004';
session_prefix = 'Session_20150403_1145';
runs           = [2:11 13:17];  % In case there are irrelevant runs recorderd to check stimulus code for presentation

%% Get EEG data
nr_runs = length(runs);   % number of runs in a session
el_data = load(fullfile(project_path,'Data', session_name, 'raw', session_prefix));
eeg_ts  = cell(1,nr_runs);

% which fields in el_data correspond to eeg data?
fields = fieldnames(el_data);
which_fields = find(~cellfun(@isempty, strfind(fields, session_prefix)));
which_fields = which_fields(runs);

% pull out eeg data from el_data structure and store in cell array eeg_ts
for ii = 1:nr_runs 
    this_field = fields{which_fields(ii)};
    eeg_ts{ii} = el_data.(this_field)'; % tranpose so that eeg_ts is time x channels
end

% pull out the impedance maps
tmp =  find(~cellfun(@isempty, strfind(fields, 'Impedances')));
impedances = cell(1, length(tmp));
for ii = 1:length(tmp), impedances{ii} = el_data.(fields{tmp(ii)}); end

clear el_data;

%% Get conditions

% Make a flicker sequence as presented in the experiment
% start_signal = eeg_make_flicker_sequence(nr_flashes, dur, isi, s_rate_eeg, 10);
load('start_signal')

% Get events file in useful units (seconds)
ev_pth = fullfile(project_path,'Data', session_name, 'raw', [session_prefix '.evt']);

% Extract the triggers from the file, and put them in timeseries
[ev_ts, start_inds] = eeg_get_triggers(ev_pth,...
    s_rate_eeg, s_rate_monitor, runs, eeg_ts, start_signal, plot_figures);

clear ev_pth start_signal;

% Find epoch onset times in samples (if we record at 1000 Hz, then also in ms)
epoch_starts = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run,...
    epochs_per_block);

% extract conditions from behavioral matfiles

directory_name = fullfile(project_path, 'Data', session_name, 'behavior_matfiles');
dir = what(directory_name);
which_mats = dir.mat(runs);

conditions = cell(1,nr_runs);
for ii = 1:nr_runs
    stimulus_file   = load(fullfile(directory_name, which_mats{ii}),'stimulus');
    sequence        = find(stimulus_file.stimulus.trigSeq > 0);
    conditions{ii}  = stimulus_file.stimulus.trigSeq(sequence)';
end

%% Epoch the EEG data
n_samples = cellfun(@length, ev_ts);
onsets    = make_epoch_ts(conditions, epoch_starts, n_samples);

epoch_time  = [0  mode(diff(epoch_starts{1}))-1]/s_rate_eeg;
ts = [];  conditions=[];
for ii = 1:nr_runs     
        [thists, this_conditions] = meg_make_epochs(eeg_ts{ii}, onsets{ii}, epoch_time, s_rate_eeg);
        ts          = cat(2,ts, thists);
        conditions  = cat(2, conditions, this_conditions);
end


%% PREPROCESS DATA
[sensorData, badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
    var_threshold, bad_channel_threshold, bad_epoch_threshold, 'eeg128xyz', verbose);

sensorData = sensorData(:,~badEpochs,~badChannels);

% Make design matrix
design = zeros(length(conditions), 3);
design(conditions==1,1) = 1; % condition 1 is full field
design(conditions==5,2) = 1; % condition 5 is right (??)
design(conditions==7,3) = 1; % condition 7 is left (??)

design     = design(~badEpochs,:);

%% Denoise data and solve GLM

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq = megGetSLandABfrequencies(0:150, .995, 12);

% denoise parameters (see denoisedata.m)
optsl.pcchoose        = -10;   % denoise with exactly 10 PCs for stimulus locked
optbb.pcchoose        = -10;   % denoise with exactly 10 PCs for broadband
optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband
optbb.npoolmethod       = {'r2','n',60};
optbb.verbose           = true;
optsl.npoolmethod       = {'r2','n',60};
optsl.verbose           = true;

% The denoise algorithm needs:
% data      : time series [channel x time samples x epoch]
% design    : design matrix [epoch x nconds]

% Permute sensorData for denoising
sensorData = permute(sensorData, [3 1 2]);

% ********* Denoise the data ********************

%   Denoise for broadband analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);

fname = fullfile(project_path, 'Data',session_name,'processed',[session_prefix '_denoisedData_bb']);

parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, ...
    'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
    'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)

%   Denoise for stimulus-locked analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);

fname = fullfile(project_path, 'Data',session_name,'processed',[session_prefix '_denoisedData_sl']);
parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, ...
    'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
    'badChannels', badChannels, 'badEpochs', badEpochs,  'opt', optsl)

%% Visualize

% Make a headplot of data before and denoising

figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(slresults.origmodel(1),icond, 'SNR');
    %clims_sl = [0, max(sl_snr1)];
    clims_sl = [0,25.6723];
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bbresults.origmodel(1),  icond, 'SNR');
    ab_snr2 = getsignalnoise(bbresults.finalmodel(1), icond, 'SNR');
    clims_ab = [0, max([ab_snr1, 12.4445])];
    %clims_ab = [0, max([ab_snr1, ab_snr2])];
    
    % convert back into 157-channel space
    ab_snr1a = to157chan(ab_snr1,~badChannels,'nans');
    ab_snr2a = to157chan(ab_snr2,~badChannels,'nans');
    sl_snr1a = to157chan(sl_snr1,~badChannels,'nans');
    
    % plot spatial maps
    subplot(3,3,(icond-1)*3+1)
    [~,ch] = megPlotMap(sl_snr1a,clims_sl,gcf,'jet',sprintf('%s : Stimulus Locked Original', condNames{icond}));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('SL no DN %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+2)
    [~,ch] = megPlotMap(ab_snr1a,clims_ab,gcf,'jet',sprintf('%s Original', condNames{icond}));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+3)
    [~,ch] = megPlotMap(ab_snr2a,clims_ab,gcf,'jet',sprintf('%s : Denoised PC %d',condNames{icond}, bbresults.pcnum(1)));
    makeprettyaxes(gca,9,9);
    makeprettyaxes(ch,9,9);
    title(sprintf('Broadband Post %s', condNames{icond}))
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure5_examplesubject3'),[],0,'.',1);
end


