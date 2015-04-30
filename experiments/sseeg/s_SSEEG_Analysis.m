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
save_data             = true;     % Save data in processed folder or not?
images_per_block      = 72;       % number of images shown within each 6-s block of experiment
epochs_per_block      = 6;        % bin data into 1 second epochs (blocks are 6 s)
blocks_per_run        = 12;       % number of blocks in one experimental run
var_threshold         = [.05 20]; % acceptable limits for variance in an epoch, relative to median of all epochs
bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
data_channels         = 1:128;    
verbose               = true;
which_subject         = 'wlsubj019';

late_timing_thresh    = 1000;     % if diff between two epoch onsets is > this value, toss the epoch 
early_timing_thresh   = 992;      % if diff between two epoch onsets is < this value, toss the epoch 


%% Define variables for this particular subject's session

switch which_subject
    case 'wlsubj019'
        session_name   = 'Session_20150417_wlsubj019';
        session_prefix = 'Session_20150417_1351';
        runs           = 1:15;  % In case there are irrelevant runs recorderd to check stimulus code for presentation
    case 'wlsubj004'
        session_name   = 'SSEEG_20150403_wl_subj004';
        session_prefix = 'Session_20150403_1145';
        runs           = [2:11 13:17]; % This is for wl_subj001: 2:9;  % In case there are irrelevant runs recorderd to check stimulus code for presentation
end
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

%% Extract conditions and initializing sequence from behavioral matfiles

directory_name = fullfile(project_path, 'Data', session_name, 'behavior_matfiles');
dir = what(directory_name);
which_mats = dir.mat(runs);

conditions  = cell(1,nr_runs);
    for ii = 1:nr_runs
          stimulus_file   = load(fullfile(directory_name, which_mats{ii}),'stimulus');
          sequence        = find(stimulus_file.stimulus.trigSeq > 0);
          conditions{ii}  = stimulus_file.stimulus.trigSeq(sequence)';      
            if isempty(strfind(session_name, '20150129'))
            else 
               conditions{ii} = conditions{ii}(1:12:end);
            end
    end

init_seq    = stimulus_file.stimulus.flashTimes.flip;
init_times  = round((init_seq - init_seq(1)) * 1000)+1;
init_ts     = ones(1, init_times(end)+40);

    for jj = 1:2:size(init_times,2)-1
            init_ts(init_times(jj):init_times(jj+1)) = false;
    end

%% Extract event times from .evt file, and define epoch onset times in samples 
%  (or in ms if you have sampled at 1000Hz)

ev_pth = fullfile(project_path,'Data', session_name, 'raw', [session_prefix '.evt']);

[ev_ts, start_inds] = eeg_get_triggers(ev_pth,...
    s_rate_eeg, s_rate_monitor, runs, eeg_ts, init_ts, plot_figures);

    if isempty(strfind(session_name, '20150129'))
        epoch_starts = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run,...
                                                                epochs_per_block);
    else
        epoch_starts = sseeg_find_epochs_pilot_eline(ev_ts, images_per_block, ...
                                                 blocks_per_run, epochs_per_block);
    end      

directory_name = fullfile(project_path, 'Data', session_name, 'behavior_matfiles');
thisdir = what(directory_name);
which_mats = thisdir.mat(runs);

clear ev_pth init_ts;    
%% Remove epochs with timing errors

    for ii = 1:nr_runs
          lag_epochs      = find(diff(epoch_starts{ii}) > late_timing_thresh)+1;
          early_epochs    = find(diff(epoch_starts{ii}) < early_timing_thresh)+1;
          time_errors     = cat(2, lag_epochs, early_epochs);
          epoch_starts{ii}(time_errors) = [];
          conditions{ii}(time_errors) = [];
    end
    
clear lag_epochs early_epochs time_errors;
%% Epoch the EEG data
n_samples = cellfun(@length, ev_ts);
onsets    = make_epoch_ts(conditions, epoch_starts, n_samples);

epoch_time  = [0  mode(diff(epoch_starts{1}))-1]/s_rate_eeg;
ts = [];  conditions=[];
    for ii = 1:nr_runs     
            [thists, this_conditions] = meg_make_epochs(eeg_ts{ii}, onsets{ii}, epoch_time, s_rate_eeg);
            ts          = cat(2, ts, thists);
            conditions  = cat(2, conditions, this_conditions);
    end

%% PREPROCESS DATA
[sensorData, badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
    var_threshold, bad_channel_threshold, bad_epoch_threshold, 'eeg128xyz', verbose);

sensorData = sensorData(:,~badEpochs,~badChannels);

%% ********* Prepare and solve GLM *********

% Make design matrix
design = zeros(length(conditions), 3);
design(conditions==1,1) = 1; % condition 1 is full field
design(conditions==5,2) = 1; % condition 5 is right (??)
design(conditions==7,3) = 1; % condition 7 is left (??)

design     = design(~badEpochs,:);

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq = megGetSLandABfrequencies((0:150)/.995, .995, 12/.995);

% denoise parameters (see denoisedata.m)
opt.pcchoose          = -10; %1.05;  % denoise with exactly 10 PCs for stimulus locked and BB
opt.npoolmethod       = {'r2','n',60};
opt.verbose           = true;
opt.pcn               = 10;
optsl = opt;
optbb = opt;
optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband


% The denoise algorithm needs:
% data      : time series [channel x time samples x epoch]
% design    : design matrix [epoch x nconds]

% Permute sensorData for denoising
sensorData = permute(sensorData, [3 1 2]);

%% ********* Denoise the data *********

%   Denoise for broadband analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);

% If requested: Save data
    if save_data
    fname = fullfile(project_path, 'Data',session_name,'processed',[session_prefix '_denoisedData']);
    parsave([fname '_bb_full70.mat'], 'results', results, 'evalout', evalout, ...
        'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
        'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)
    end

  Denoise for stimulus-locked analysis
[results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);

    if save_data
    parsave([fname '_sl_full.mat'], 'results', results, 'evalout', evalout, ...
        'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
        'badChannels', badChannels, 'badEpochs', badEpochs,  'opt', optsl)
    end

data_orig = zeros(1,128);
data_final = zeros(1,128);
tmp = find(badChannels);
data_orig(~badChannels) = results.origmodel.r2;
data_final(~badChannels) = results.finalmodel.r2;
figure(1); plotOnEgi(data_orig); title('Original'); colorbar; 
figure(2); plotOnEgi(data_final); title('Final'); colorbar; 
figure(3); plotOnEgi(data_final - data_orig); title('Final minus Original'); colorbar; 

%% Visualize
% sseegMakePrePostHeadplot(project_path,session_name,session_prefix,true)
