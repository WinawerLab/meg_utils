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
% nr_flashes = 4; 
% dur = 1/s_rate_monitor; 
% isi = 7/s_rate_monitor; 

%% Define variables for this subject's session
session_name   = 'SSEEG_20150403_wl_subj004';
session_prefix = 'Session_20150403_1145';
runs           = [2:11 13:17]; % In case there are irrelevant runs recorderd to check stimulus code for presentation

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
% start_signal = eeg_make_flicker_sequence(nr_flashes, dur, isi, s_rate_eeg, 10);
load('start_signal')

% Get events file in useful units (seconds)
ev_pth = fullfile(project_path,'Data', session_name, 'raw', [session_prefix '.evt']);

% Extract the triggers from the file, and put them in timeseries
[ev_ts, start_inds, t] = eeg_get_triggers(ev_pth,...
    s_rate_eeg, s_rate_monitor, runs, eeg_ts, start_signal, plot_figures);

%% Find epochs

samples_per_epoch = 1000;
epoch_starts = sseeg_find_epochs(ev_ts, trigs_per_block, blocks_per_run, DINs_per_epoch);

%% extract conditions from behavioral matfiles

directory_name = fullfile(project_path, 'Data', session_name, 'behavior_matfiles');
dir = what(directory_name);
which_mats = dir.mat(runs);

order_long = cell(1,nr_runs);
for ii = 1:nr_runs
    stimulus_file   = load(fullfile(directory_name, which_mats{ii}),'stimulus');
    sequence        = find(stimulus_file.stimulus.trigSeq > 0);
    order_long{ii}  = stimulus_file.stimulus.trigSeq(sequence)';
end

%% create inputs necessary for eeg_make_epochs function

% order       = [1 3 1 3 5 3 5 3 7 3 7 3]; % this parameter might go to the top of script
epoch_ts    = make_epoch_ts(order_long, nr_runs, ev_ts, epoch_starts);

%% run eeg_make_epochs

ts_cell         = cell(1,nr_runs);
conditions_all  = cell(1,nr_runs);
for ii = 1:nr_runs
    epoch_time  = [1  mode(diff(epoch_starts{1}))+1]; 
        [ts, conditions] = eeg_make_epochs(eeg_ts{ii}', epoch_ts{ii}, epoch_time, s_rate_eeg);
             ts_cell{ii} = ts;
      conditions_all{ii} = conditions;
end

%% ***** just playing around with the ecogCalcOnOffSpectra function ******
%  *****               will be deleted eventually                   ******

% off.signal = ts_cell{ii}(:, find(conditions_all{ii} == 3),83);
% off.signal = off.signal(:,1:24);
% full  = find(conditions_all{ii} == 1);
% right = find(conditions_all{ii} == 5);
% left  = find(conditions_all{ii} == 7);
% on.signal = ts_cell{ii}(:, [full left], 83);
% 
% [on, off] = ecogCalcOnOffSpectra(on, off, 1, 0);
% 
% t = [.001:.001:.995];
% stimF = 12;
% calcPower = 1;
% fH = ecogPlotOnOffSpectra(on, off, t, stimF, calcPower)
% ave_on = mean(on.meanFFT,3);
% ave_off = mean(off.meanFFT,3);
% 
% figure; plot(loglog(ave_on)); 
% hold on; 
% plot(loglog(ave_off));
% axis([5 50 0 3]);
