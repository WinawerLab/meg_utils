%% Script to analyze EEG pilot data

% Description comes here: Script to analyze a pilot steady state EEG data
% set (subject wl_s004). Contrast patterns were contrast-reversed at 12 Hz
% in 6 second blocks alaternating with 6 s of blank (mean luminance) while
% subjects fixated the middle of the screen and detected a fixation color
% change. Stimuli consisted of either full-field (11 deg radius?? check
% this), left field, or right field apertures.


% Dependencies: meg_utils github repository


%% Define variables for this experiment
project_path        = '/Volumes/server/Projects/EEG/SSEEG/';
s_rate_eeg          = 1000; % sample rate of the eeg in Hz
s_rate_monitor      = 60;   % sample rate of the monitor in Hz
plot_figures        = true; % Plot debug figures or not?
images_per_block    = 72;   % number of images shown within each 6-s block of experiment
triggers_per_epoch  = 6;    % black-to-white transitions
blocks_per_run      = 6;    % number of blocks in one experimental run
remove_bad_epochs   = true;

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

%% Find epoch onset times in samples (if we record at 1000 Hz, then also in ms)
epoch_starts = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run);

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
<<<<<<< HEAD
epoch_ts    = make_epoch_ts(conditions, nr_runs, ev_ts, epoch_starts);
=======

% order       = [1 3 1 3 5 3 5 3 7 3 7 3]; % this parameter might go to the top of script

epoch_ts    = make_epoch_ts(order_long, nr_runs, ev_ts, epoch_starts);
>>>>>>> 9ea37560b20416cbcdabc25b5003f85e6b29d53e

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

%% CALCULATE FOURIER TRANSFORMS

t                   = size(ts,1);
num_epoch_time_pts  = size(ts,1);

full                = find(conditions == 1);
right               = find(conditions == 5);
left                = find(conditions == 7);
off                 = find(conditions == 3);

ts_off_full         = ts(:, off, data_channels);
ts_on_full          = ts(:, full, data_channels);
ts_on_right         = ts(:, right, data_channels);
ts_on_left          = ts(:, left, data_channels);
freq                = (0:num_epoch_time_pts-1)/(num_epoch_time_pts/1000);
channelsToPlot      = 90; 

ft_on_epoched_full   = fft(ts_on_full) / length(t)*2;   
ft_on_epoched_right  = fft(ts_on_right) / length(t)*2;
ft_on_epoched_left   = fft(ts_on_left) / length(t)*2; 
ft_off_epoched_full  = fft(ts_off_full)/ length(t)*2; 

amps_on_full  =  abs(ft_on_epoched_full);   clear ft_on_epoched_full;
amps_on_right  =  abs(ft_on_epoched_right);  clear ft_on_epoched_right;
amps_on_left  =  abs(ft_on_epoched_left);   clear ft_on_epoched_left;
amps_off_full = abs(ft_off_epoched_full);   clear ft_off_epoched_full;

    on_full = squeeze(nanmedian(amps_on_full(:,:,channelsToPlot), 2));
    on_right = squeeze(nanmedian(amps_on_right(:,:,channelsToPlot), 2));
    on_left = squeeze(nanmedian(amps_on_left(:,:,channelsToPlot), 2));
    off_full = squeeze(nanmedian(amps_off_full(:,:,channelsToPlot), 2));

    mean_on_full = mean(on_full,2);
    mean_on_right = mean(on_right,2);
    mean_on_left = mean(on_left,2);
    mean_off_full = mean(off_full,2);
    
  figure(104);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20)
    hold on;
    plot(freq, off_full, 'LineWidth', 2);
    xlim([3 85])
    yl = get(gca, 'YLim');
    % for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (Tesla)');
    title('Stimulus On');
    get_MEG_axes('True');
    
        
  figure(205);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20)
    hold on;
    plot(freq, mean_off, 'LineWidth', 2)
    xlim([3 85])
    yl = get(gca, 'YLim');
    % for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (Tesla)')
    title('Blank Periods')
    get_MEG_axes('True');


[amps_on_full,amps_on_right,amps_on_left, ...
    amps_off_full,amps_off_right,amps_off_left] = ssmeg_fourier(...
    t, num_epoch_time_pts, ts_on_full, ts_off_full, num_epochs, [], [], 1)

%% Alternative method to extract conditions 

<<<<<<< HEAD
el_data         = load(fullfile(project_path,'Data', session_name, 'raw', session_prefix));
condition_cell  = el_data.ECI_TCPIP_55513;
conditions      = str2double(condition_cell(1,:));
clear el_data;


=======
% off.signal = ts_cell{ii}(:, find(conditions == 3), 81:85);
% off.signal = off.signal(:,1:24,:);
% full  = find(conditions == 1);
% right = find(conditions == 5);
% left  = find(conditions == 7);
% on.signal = ts_cell{ii}(:, [full left], 81:85);
% 
% [on, off] = ecogCalcOnOffSpectra(on, off, 1, 0);
% ave_on = mean(on.meanFFT,3);
% ave_off = mean(off.meanFFT,3);
% 
% figure; plot(loglog(ave_on)); 
% hold on; 
% plot(loglog(ave_off));
% axis([5 50 0 3]);
>>>>>>> 9ea37560b20416cbcdabc25b5003f85e6b29d53e

