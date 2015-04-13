%% Script to analyze EEG pilot data

% Description comes here: Script to analyze a pilot steady state EEG data
% set (subject wl_s004). Contrast patterns were contrast-reversed at 12 Hz
% in 6 second blocks alaternating with 6 s of blank (mean luminance) while
% subjects fixated the middle of the screen and detected a fixation color
% change. Stimuli consisted of either full-field (11 deg radius?? check
% this), left field, or right field apertures.


% Dependencies: meg_utils github repository


%% Define variables for this experiment
project_path          = '/Volumes/server/Projects/EEG/SSEEG/';
s_rate_eeg            = 1000; % sample rate of the eeg in Hz
s_rate_monitor        = 60;   % sample rate of the monitor in Hz
plot_figures          = true; % Plot debug figures or not?
images_per_condition  = 72;   % number of images shown within each 6-s block of experiment
epochs_per_condition  = 6;    % black-to-white transitions
conditions_per_run    = 12;   % number of blocks in one experimental run
remove_bad_epochs     = true; 

% Photodiode start sequence parameters
%   We look for this in the trigger sequence to indicate where the
%   experiment starts. In this case, we flash the square for one frame
%   every 8 frames, 4x
%  nr_flashes = 4; 
%  dur = 1/s_rate_monitor; 
%  isi = 7/s_rate_monitor; 

%% Define variables for this subject's session
session_name   = 'SSEEG_20150403_wl_subj004';
session_prefix = 'Session_20150403_1145';
runs           = [2:11 13:17];  % In case there are irrelevant runs recorderd to check stimulus code for presentation

%% Get toolboxes and code
% addpath(fullfile(project_path, 'Code'));
addpath(genpath('~/matlab/git/meg_utils'));

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
    eeg_ts{ii} = el_data.(this_field);
end

% pull out the impedance maps
tmp =  find(~cellfun(@isempty, strfind(fields, 'Impedances')));
impedances = cell(1, length(tmp));
for ii = 1:length(tmp), impedances{ii} = el_data.(fields{tmp(ii)}); end

clear el_data;

%% Get timeseries from event files

% Make a flicker sequence as presented in the experiment
% start_signal = eeg_make_flicker_sequence(nr_flashes, dur, isi, s_rate_eeg, 10);
load('start_signal')

% Get events file in useful units (seconds)
ev_pth = fullfile(project_path,'Data', session_name, 'raw', [session_prefix '.evt']);

% Extract the triggers from the file, and put them in timeseries
[ev_ts, start_inds] = eeg_get_triggers(ev_pth,...
    s_rate_eeg, s_rate_monitor, runs, eeg_ts, start_signal, plot_figures);

clear ev_pth start_signal;
%% Find epoch onset times in samples (if we record at 1000 Hz, then also in ms)
epoch_starts = sseeg_find_epochs(ev_ts, images_per_condition, conditions_per_run,...
    epochs_per_condition);

%% extract conditions from behavioral matfiles

directory_name = fullfile(project_path, 'Data', session_name, 'behavior_matfiles');
dir = what(directory_name);
which_mats = dir.mat(runs);

conditions = cell(1,nr_runs);
for ii = 1:nr_runs
    stimulus_file   = load(fullfile(directory_name, which_mats{ii}),'stimulus');
    sequence        = find(stimulus_file.stimulus.trigSeq > 0);
    conditions{ii}  = stimulus_file.stimulus.trigSeq(sequence)';
end

%% create inputs necessary for eeg_make_epochs function

epoch_ts    = make_epoch_ts(conditions, nr_runs, ev_ts, epoch_starts);

%% run eeg_make_epochs

epoch_time  = [1  mode(diff(epoch_starts{1}))+1];
ts = [];  conditions=[];
for ii = 1:nr_runs     
        [thists, this_conditions] = eeg_make_epochs(eeg_ts{ii}', epoch_ts{ii}, epoch_time, s_rate_eeg);
        ts          = cat(2,ts, thists);
        conditions  = cat(2, conditions, this_conditions);
end

%% PREPROCESS DATA

data_channels = 1:128;

if remove_bad_epochs
    
    % This identifies any epochs whos variance is outside some multiple of the
    % grand variance
    bad_epochs = meg_find_bad_epochs(ts(:,:,data_channels), [.05 20]);
    
    % any epoch in which more than 10% of channels were bad should be removed
    % entirely
    epochs_to_remove = mean(bad_epochs,2)>.2;
    
    % once we remove 'epochs_to_remove', check whether any channels have more
    % than 10% bad epochs, and we will remove these
    channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.2;
    
    bad_epochs(epochs_to_remove,:) = 1;
    bad_epochs(:,channels_to_remove) = 1;
    
    figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')
    
    rawts =ts;
    ts = meg_remove_bad_epochs(bad_epochs, rawts);
end

%% Calculate absolute values of Fourier transformed data

t                   = size(ts,1);
num_epoch_time_pts  = size(ts,1);
channels_to_plot    = 50:99;
num_epochs          = epochs_per_condition;

[amps_on_full,amps_on_right,amps_on_left, amps_off_full,amps_off_right, ...
    amps_off_left] = sseeg_fourier(t, num_epoch_time_pts, ts(:,:,data_channels), conditions, ...
        num_epochs, data_channels, channels_to_plot, 1); 
    
%% Plot values on mesh

[coh_full, coh_right, coh_left] = sseeg_cal_coh(amps_on_full, amps_on_right, ...
     amps_on_left, [], 1, 1);
    
%% Calculate broadband spectra **** UNFINISHED

off_conditions      = find(conditions ==3);
a                   = size(off_conditions,2)/3;

off_full  = ts(:,off_conditions(1:a), :);
off_right = ts(:,off_conditions(a+1:2*a), :);
off_left  = ts(:,off_conditions((2*a)+1:3*a), :);
on_full   = ts(:, find(conditions == 1), :);   
on_right  = ts(:, find(conditions == 5), :);
on_left   = ts(:, find(conditions == 7), :);   

amps = zeros(size(amps_on_full,1), size(amps_on_full,2), size(amps_on_full,3), 6); 
amps(:,:,:,1) = on_full;
amps(:,:,:,2) = on_right;
amps(:,:,:,3) = on_left;
amps(:,:,:,4) = off_full;
amps(:,:,:,5) = off_right;
amps(:,:,:,6) = off_left;

[ab, log_frequency, bbPoly, frequencies] = sseeg_broadband(num_epoch_time_pts, amps);

%% Alternative method to extract conditions 

% el_data         = load(fullfile(project_path,'Data', session_name, 'raw', session_prefix));
% condition_cell  = el_data.ECI_TCPIP_55513;
% conditions      = str2double(condition_cell(1,:));
% clear el_data;

