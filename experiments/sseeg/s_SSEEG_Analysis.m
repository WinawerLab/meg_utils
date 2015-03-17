%% Script to analyze EEG data

% Description comes here: ...


% Note: Vistadisp code starts with a white square in the upper left corner
% to trigger photodiode.
%
% Dependencies: meg_utils github repository


%% Define variables for this experiment
project_path     = '/Volumes/server/Projects/EEG/SSEEG/';
s_rate_eeg       = 1000; % sample rate of the eeg in Hz
s_rate_monitor   = 60;   % sample rate of the monitor in Hz
num_time_pts     = 1000; % Number of time points for 1 epoch
plot_figures     = true; % Plot debug figures or not?
trigs_per_block  = 72;   % number of contrast reversals in one block of experiment
blocks_per_run   = 6;    % number of blocks in one experimental run
DINs_per_epoch   = 6;
% Photodiode start sequence parameters
%   We look for this in the trigger sequence to indicate where the
%   experiment starts. In this case, we flash the square for one frame
%   every 8 frames, 4x
nr_flashes = 4; 
dur = 1/s_rate_monitor; 
isi = 7/s_rate_monitor; 

%% Define variables for this subject's session
session_name   = 'Pilot_SSEEG_20150129_wl_subj001';
session_prefix = 'Session_20150129_1007';
runs             = 3:10; % In case there are irrelevant runs recorderd to check stimulus code for presentation

%% Get toolboxes and code
addpath(fullfile(project_path, 'Code'));
addpath(genpath('~/matlab/git/meg_utils'));

%% Get EEG data
nr_runs = length(runs);   % number of runs in a session
el_data = load(fullfile(project_path,'Data', session_name, 'raw', session_prefix));
eeg_ts  = cell(1,nr_runs);

% which fields in el_data correspond to eeg data?
fields = fieldnames(el_data);
which_fields = find(~cellfun(@isempty, strfind(fields, session_prefix)));
which_fields = which_fields(runs);

for ii = 1:nr_runs % Assuming that there will be a DIN, TCPIP and Marks field
    this_field = fields{which_fields(ii)};
    eeg_ts{ii} = el_data.(this_field);
end

clear el_data;

%% Get timeseries from event files

% Make a flicker sequence as presented in the experiment
start_signal = eeg_make_flicker_sequence(nr_flashes, dur, isi, s_rate_eeg, 10);

% Get events file in useful units (seconds)
ev_pth = fullfile(project_path,'Data', session_name, 'raw', [session_prefix '.evt']);

% Extract the triggers from the file, and put them in timeseries
[ev_ts, start_inds, t] = eeg_get_triggers(ev_pth,...
    s_rate_eeg, s_rate_monitor, runs, eeg_ts, start_signal, plot_figures);

%% Find epochs

samples_per_epoch = 1000;
epoch_starts = sseeg_find_epochs(ev_ts, trigs_per_block, blocks_per_run, DINs_per_epoch);

%% prepare inputs for meg_make_epochs function

triggers = cell(1,nr_runs);
for ii = 1:nr_runs
     triggers{ii} = zeros(1,length(ev_ts{ii}));
     triggers{ii}(round(on_epochs{ii}(1:12)*s_rate_eeg))  = 1;
     triggers{ii}(round(on_epochs{ii}(13:24)*s_rate_eeg)) = 3;
     triggers{ii}(round(on_epochs{ii}(25:36)*s_rate_eeg)) = 5;
     triggers{ii}(round(off_epochs{ii}(1:36)*s_rate_eeg)) = 2;
end

epoch_time = [on_epochs{ii}(1) on_epochs{ii}(2)]; 
fs = s_rate_eeg;
raw_ts = eeg_ts{ii}';
trigger = triggers{ii};

[ts, conditions] = meg_make_epochs(raw_ts, trigger, epoch_time, s_rate_eeg);
%% Alternative for loop in the meg_make_epochs function

for ii = 1:num_epochs
    inds = onsets(ii):onsets(ii)+(epoch_len-1);
    ts(ii, :, :) = raw_ts(inds,:);    
end

