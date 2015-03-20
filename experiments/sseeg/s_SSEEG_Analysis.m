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

%% create inputs necessary for eeg_make_epochs function

order       = [1 3 1 3 5 3 5 3 7 3 7 3]; % this parameter might go to the top of script
epoch_ts    = make_epoch_ts(order, nr_runs, ev_ts, epoch_starts);

%% run eeg_make_epochs
ts_cell = cell(1,nr_runs);

for ii = 1:nr_runs
    epoch_time  = [epoch_starts{ii}(1) epoch_starts{ii}(2)]; 
        [ts, conditions] = eeg_make_epochs(eeg_ts{ii}', epoch_ts{ii}, epoch_time, s_rate_eeg);
    ts_cell{ii} = ts;
end

%%
off.signal = ts(:, find(conditions == 3), 70);

full  = find(conditions == 1);
right = find(conditions == 5);
left  = find(conditions == 7);
on.signal = ts(:, [full right], 70);


[on, off] = ecogCalcOnOffSpectra(on, off, 1, 1)
