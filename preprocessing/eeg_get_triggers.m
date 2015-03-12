function [ev_ts, t] = eeg_get_triggers(ev_pth, s_rate_eeg, ...
    s_rate_monitor, runs, eeg_ts, start_signal, plot_figures)

%% Function to extract triggers from photodiode responses saved in evt file

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
    num_time_pts = size(eeg_ts{ii},2);
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
    start_ind = eeg_find_trial_begin(start_signal, ev_ts{which_run});
    ev_ts{which_run}(1:start_ind) = 0;

    if plot_figures, start_inds(which_run) = start_ind; end %#ok<AGROW>
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

%% Find on trigger times
% Find all of the reversals from 0-->1 or 1-->0 in ev_ts. Then, from the
% start index of the run, we make 12th reversal an epoch timepoint, since
% there are two reversals per DIN, and six DIN's per epoch. A couple hacks,
% start_run_ind{ii} is shifted by one from reversal_inds. Additionally, we
% might want a better way of defining the value '431', it is 432-1, 432
% being 72 times 6 (reversals per block times blocks
% per run).

trigs_per_run = trigs_per_block * blocks_per_run;
reversal_inds  = cell(1,nr_runs);
on_epoch_times = cell(1,nr_runs);

for ii = 1:nr_runs
    reversal_inds{ii} = find(diff(ev_ts{ii}==1));
    start_inds = find(reversal_inds{ii}+1 == start_run_ind{ii});
    on_epoch_times{ii} = t{ii}(reversal_inds{ii}(start_inds:12:start_inds+(trigs_per_run*2)-1));
end

if plot_figures
       hold on; plot(on_epoch_times{ii}, 1, 'rx');
end

return
%% Find off trigger times, not working yet

off_epoch_times = cell(1,nr_runs);
    big_diffs = find(diff(reversal_inds{ii}) > 1000);
    off_inds = big_diffs(big_diffs > start_inds(1,ii));
for ii = first_run:nr_runs
    
    
    
%         for ll = 1:6
%         off_period_inds( = reversal_inds{ii}(off_inds(ll)):reversal_inds{ii}(off_inds(ll))+5000
%         off_epoch_times{ii} = t{ii}(off_period_inds(ll:1000:ll+))
end


%%

% % Transform indices into trigger times.
% for which_run = 1:nr_runs
%     start_run_times{which_run} = t{which_run}(start_run_ind{which_run});  
% end
% 
% % Take only every 6th DIN trigger per run
% for ii = 1:nr_runs
%     triggers{ii} = ev_ts{ii}(start_run_times{ii}:6:end)
% end
% 
% % Return triggers, ev_ts and t to main script.
% 
% 
% 
% return


%% STILL WORK IN PROGRESS FROM HERE, NOT SURE WHAT THIS IS

% 
% triggers_start_time = [];
% triggers_start_ind = [];
% for ii = 1:size(start_run_time,2)
%     temp = find(t(ev_inds) == start_run_time(ii));
%     
%     triggers_start_ind{ii} = ev_inds(temp+1);
%     triggers_start_time{ii} = t(ev_inds(temp+1));
%     
% end
% 
% triggers_start_ind = cell2mat(triggers_start_ind);
% triggers_start_time = cell2mat(triggers_start_time);
% triggers_start_time = triggers_start_time(first_run:end);
% triggers_start_ind = triggers_start_ind(first_run:end);
% % clear start_run_ind    start_run_time    end_init_seq    init_seq_diff_ind; 
% 
% %% Plot
%      if plot_figures
%          % Plot only relevant data
%          figure; plot(t, ev_ts); hold on; plot(t,eeg_ts(1:4,:));
%          plot(triggers_start_time, 100, 'r+');
%          xlabel('Time [ms]')
%          ylabel('Electric field strength (microVolts and arbitrary units for triggers'); 
%      end
% 
% %% From the triggers, work forwards, designating relevant epochs and on/off periods 
% % nr_cycles_in_epoch = 6; %ms
% % 
% % triggers = ev_ts(triggers_start_ind(1):triggers_start_ind(1));
% % 
% % 
% % bins = triggers_start_ind(1):1000:triggers_start_ind(2)-1;
% % figure; plot(t, ev_ts); hold on; plot(t(bins),100*ones(length(t(bins)),1), 'ro');
% % 
% % times = t(higher_than_rf);
% % 
% % 
% % bins_2 = triggers_start_ind(1):1000:triggers_start_ind(2)-1;
% % 
% % 
% % figure; plot(t, ev_ts); hold on;  plot(t(higher_than_rf), 100*ones(length(t(higher_than_rf)),1), 'ro');
% 
% %% Define on epochs as every six triggers
% epoch_size = 1;
% blocks_per_run = 6;
% flips_per_block = 72;
% reversals_per_block = flips_per_block*blocks_per_run;
% nr_runs = size(triggers_start_time,2);
% 
%         periods_on_ind = diff(ev_ts(triggers_start_time(1):triggers_start_time(2)-1)==100);
%         frame_times    = sort([find(periods_on_ind==1)';find(periods_on_ind==-1)']);
%     
% on_epochs = zeros(36, nr_runs);
% for ii = 1:nr_runs
%     
%     if ii < nr_runs
%         periods_on_ind = diff(ev_ts(triggers_start_time(ii):triggers_start_time(ii+1)-1)==100);
%     else 
%         periods_on_ind = diff(ev_ts(triggers_start_time(ii):t(end))==100);
%     end
%     frame_times = sort([find(periods_on_ind==1)';find(periods_on_ind==-1)']);
%     on_epochs(:,ii) = t(frame_times(1:12:blocks_per_run*flips_per_block))+triggers_start_time(ii);    
% 
% end
% 
% %% Find off epochs
% 
% off_epochs = zeros(6, nr_runs);
% for ii = 1:nr_runs
% 
%     if ii < nr_runs
%         periods_on_ind = diff(ev_ts(triggers_start_time(ii):triggers_start_time(ii+1)-1)==100);
%     else 
%         periods_on_ind = diff(ev_ts(triggers_start_time(ii):t(end))==100);
%     end
%         frame_times = sort([find(periods_on_ind==1)';find(periods_on_ind==-1)']);
%         off_start_inds = find(diff(frame_times) > flips_per_block*2);
%         off_epochs(:,ii) = t(frame_times(off_start_inds(1:6))) + triggers_start_time(ii);
% 
% end 
%             
% 
% off_start_inds = find(diff(frame_times) > flips_per_block*2);
% off_start_times = t(frame_times(off_start_inds(1:6))) + triggers_start_time(1)
% size_off_epoch = (on_epochs(1*epochs_per_block+1) - off_start_times(1))/epochs_per_block;
% 
% %% plot
%      if plot_figures
%               hold on;
%               for ii = 1:size(on_epochs,2)
%                     plot(on_epochs(:,ii), 100, 'rx');
%                     plot(off_epochs(:,ii), 99, 'kx');
%               end 
%      end
