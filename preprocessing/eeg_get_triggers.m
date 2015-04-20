function [ev_ts, start_inds, t] = eeg_get_triggers(ev_pth, s_rate_eeg, ...
    s_rate_monitor, runs, eeg_ts, init, plot_figures)
%% Extract triggers from photodiode responses saved in NetStation evt file
%
% [ev_ts, t] = eeg_get_triggers(ev_pth, s_rate_eeg, ...
%       s_rate_monitor, first_run, flicker_sequence, plot_figures);
%
% Description: ...

% INPUTS:
% ev_pth            : String defining the path where you can find the .evt file
% s_rate_eeg        : Sample rate of the EEG recording (Hz)
% s_rate_monitor    : Refresh rate of the monitor used during EEG experiment
%                       in Hz
% runs              : In case there are dummy runs or practice runs, define 
%                       which runs to extract triggers from. Length of runs
%                       should match length of eeg_ts, as there must be a
%                       corresponding eeg time series for each run with
%                       triggers
% eeg_ts            : cell array of eeg data, where each cell is time x channels
%         
% start_signal      : A vector indicating the flicker sequence which marks the start of each run
% plot_figures      : Plot debug figures or not?

% OUTPUTS:
% ev_ts             : Timeseries (vector of zeros and ones) with the same 
%                       length as the eeg data per run
% t                 : Cell array where each cell is the time vector per run
%                       in seconds (same length as ev_ts and eeg_ts)

% Example: ....

% Dependencies: ...

%% Read in the event file 

% TO DO: Use different read function, for example 'fread'.
ev_times = eeg_read_evt_file(ev_pth); % event times in seconds of each DIN

%% Find start and end points of runs. 
start_indices = (find((diff(ev_times) < 0))+1)';
start_indices = [1 start_indices]; % no negative diff for the first run%% Make triggers increase continuously in time by accounting for zero reset
end_indices   = [start_indices(2:end)-1 length(ev_times)];

start_indices = start_indices(runs);
end_indices   = end_indices(runs);
nr_runs       = length(start_indices);

% chop up ev_times by run
ev_times_by_run = cell(nr_runs,1);
for ii = 1:nr_runs
     ev_times_by_run{ii} = ev_times(start_indices(ii):end_indices(ii));
end


%% Convert an array of time points to a time series.
%   The time series contains the values 1 and 0, for white square and
%   black square respectively. 

ev_ts = cell(1, nr_runs);
t     = cell(1, nr_runs);
for ii = 1:nr_runs
    num_time_pts = size(eeg_ts{ii},1);
    t{ii} = (0:num_time_pts-1) / s_rate_eeg; 
    
    [~, inds] = findnearest(ev_times_by_run{ii}, t{ii});
    ev_ts{ii} = false(size(t{ii}));
    
    
    % pad the ev_ts so that each din lasts as long as one screen refresh
    for jj = 0:ceil(s_rate_eeg/s_rate_monitor)-1
        ev_ts{ii}(inds+jj) = true;
    end
    % crop ev_ts so it is the same length as eeg_ts
    ev_ts{ii} = ev_ts{ii}(1:num_time_pts);
end
   

%% Locate initializing sequence, and clip triggers prior to experiment start
for which_run = 1:nr_runs
    start_ind = eeg_find_trial_begin(init, ev_ts{which_run});
    ev_ts{which_run}(1:start_ind) = 0;

    if plot_figures, start_inds(which_run) = start_ind; else start_inds(which_run) = start_ind; end %#ok<AGROW>
end

%% If requested, plot figures
if plot_figures
    % Plot events in all runs aligned to start time within each run
    figure;  hold all;
    for ii = 1:nr_runs
        these_inds = start_inds(ii):length(t{ii});
        plot(t{ii}(these_inds)-t{ii}(these_inds(1)), ev_ts{ii}(these_inds)+ii, 'LineWidth', 2);
    end      
end

return


end