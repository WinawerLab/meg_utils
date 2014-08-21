function [t, ts, num_epoch_time_pts,epochs_on,epochs_off,epochs_condition] = ...
    ssmeg_fix_triggers_v2(ts, experiment_name, produce_figures)

% This function will fix the misalignment of the triggers. 
%
% The trigger numbers can be build with the binary number system. In case
% of the SSMEG experiment we used 8 triggers numbers, which could be defined by the
% combination of activity in the 3 trigger channels in pairs of two.
% 
% In this case, weird triggers will have values 4 and 6.

% This function will:
% a. Plots the triggers
% b. Fixes the triggers for Full-Left-Right SSMEG experiments
% c. Prints the amount of triggers per trigger-number
% d. Defines the start of single epochs per condition and saves this order

% INPUTS:
% ts:           timeseries
% stim_file:    Parameters of the experiment in Matlab run-file

% OUTPUTS:
% num_epoch_time_pts:   Number of timepoints per single epoch (i.e. 167)
% epoch_start_*_on/off_ind: Index of start particular epoch in timeseries 


%%

% Get sampling rate
Fs = 1000; % Hz

% Get 'real' refresh rate from stimulus parameters file
% frame_dur = nanmedian(diff(stim_file.response.flip))/5*1000; % Not used anymore!

%Channels of interest
trigger_channel    = [161 162 163 164];

% Not used right know, but if needed:
%photodiode_channel = 192;
%noise_channels     = [158 159 160];

%% Get peaks for every trigger
t = (0:size(ts,1)-1) / Fs; % seconds

num_runs = 6;


% onsets trigger 1
data_trig1   = diff(ts(:,trigger_channel(1)));
thresh1      = (max(data_trig1) - min(data_trig1)) / 50;
[trig1.value, trig1.ind] = findpeaks([0; -data_trig1],'MINPEAKHEIGHT',thresh1, 'MINPEAKDISTANCE', 10);
peaks_per_run1 = length(trig1.ind) / num_runs;
start_peaks_inds1 = trig1.ind(1:peaks_per_run1:length(trig1.ind));


% onset trigger 2
data_trig2   = diff(ts(:,trigger_channel(2)));
thresh2      = (max(data_trig2) - min(data_trig2)) / 50;
[trig2.value, trig2.ind] = findpeaks([0; -data_trig2],'MINPEAKHEIGHT',thresh2, 'MINPEAKDISTANCE', 10);

peaks_per_run2 = length(trig2.ind) / num_runs;
start_peaks_inds2 = trig2.ind(1:peaks_per_run2:length(trig2.ind));


% onset trigger 3
data_trig3   = diff(ts(:,trigger_channel(3)));
thresh3      = (max(data_trig3) - min(data_trig3)) / 50;
[trig3.value, trig3.ind] = findpeaks([0; -data_trig3],'MINPEAKHEIGHT',thresh3, 'MINPEAKDISTANCE', 10);

peaks_per_run3 = length(trig3.ind) / num_runs;
start_peaks_inds3 = trig3.ind(1:peaks_per_run3:length(trig3.ind));

% % onset trigger 4
% data_trig4   = diff(ts(:,trigger_channel(4)));
% thresh4      = (max(data_trig4) - min(data_trig4)) / 50;
% [trig4.value, trig4.ind] = findpeaks([0; -data_trig4],'MINPEAKHEIGHT',thresh4, 'MINPEAKDISTANCE', 10);
% 
% peaks_per_run4 = length(trig4.ind) / num_runs;
% start_peaks_inds4 = trig4.ind(1:peaks_per_run4:length(trig4.ind));


triggerBinary = zeros(length(t), 3);
triggerBinary(trig1.ind,1) = 1;
triggerBinary(trig2.ind,2) = 1;
triggerBinary(trig3.ind,3) = 1;
% triggerBinary(trig4.ind,4) = 1;

trigger = triggerBinary * [1 2 4]';

clear triggerBinary; 
clear data_trig1; clear data_trig2; clear data_trig3; clear data_trig4;


%% Find the weird peaks
if produce_figures

    figure(200); clf
    set(gca, 'ColorOrder', jet(4), 'Color', [.6 .6 .6]); hold all; 
    plot(bsxfun(@plus, ts(:,161:164), (1:4)*40000), '-o',  'LineWidth', 2); 
    legend(cellstr(num2str((161:164)')))
% 
    weird = find(trigger == 9 | trigger == 10 | trigger == 11 | trigger == 12 | ...
    trigger == 13 | trigger == 14 | trigger == 15);
    set(gca, 'XTick', weird, 'XGrid', 'on');
    title('Trigger channels (indices and values)');
    xlabel('Time [s]');
    ylabel('Relative amplitude');
    legend('Trigger 1 / Chan 161', 'Trigger 2 / Chan 162', 'Trigger 3 / Chan 163' , 'Trigger 4 / Chan 164');
    
else
    weird = find(trigger == 9 | trigger == 10 | trigger == 11 | trigger == 12 | ...
        trigger == 13 | trigger == 14 | trigger == 15);
end



%% Fix triggers when they have value 9 or higher 
% 
for ii = 1:length(weird)
    
    % Get index of three triggers
    a = find(trig1.ind == weird(ii));
    b = find(trig2.ind == weird(ii));
    c = find(trig3.ind == weird(ii));
    c = find(trig4.ind == weird(ii));
    
    % If there is no element for trigger 1 and trigger 2, it means that
    % trigger 3 is 1 ms early. So we add +1
    if isempty(a) && isempty(b)
        trig3.ind(c) = trig3.ind(c) + 1;
    
    % If trigger 2 and trigger 3 have contain both elements, it means that
    % trigger 1 is 1 ms late. So we discard 1 ms of trigger 1.
    elseif ~isempty(b) && ~isempty(c)
        trig1.ind(find(trig1.ind==weird(ii)+1)) = trig1.ind(find(trig1.ind==weird(ii)+1)) - 1;
        
    end
end
   

   
%% Recompute triggers to check weird peaks again

triggerBinary = zeros(length(t), 3);
triggerBinary(trig1.ind,1) = 1;
triggerBinary(trig2.ind,2) = 1;
triggerBinary(trig3.ind,3) = 1;
% triggerBinary(trig4.ind,4) = 1;

trigger = triggerBinary * [1 2 4]'; clear triggerBinary;

%% check that we have the right number of triggers
fprintf('%s\n','Look at sum of trigger to check misalignment')
fprintf('%s\t\t%s\n', 'Trigger nr', 'Sum')
for ii = 1:15; fprintf('%d\t%d\n', ii, sum(trigger == ii)); end

%% Fix other weird triggers (i.e. too much of a certain trigger)

trigger_times = find(trigger);
trigger_intervals = diff(trigger_times);
[sorted, inds] = sort(trigger_intervals);

weird_triggers = find(sorted==1);
nr_weird = length(weird_triggers);

weird2 = zeros(nr_weird,1);

for ii = 1:nr_weird
    weird2(ii) = trigger_times(inds(ii));
end

for ii = 1:length(weird2);
    a = find(trig1.ind == weird2(ii));
    b = find(trig2.ind == weird2(ii));
    c = find(trig3.ind == weird2(ii));
%     d = find(trig4.ind == weird2(ii));

    % If there is no element for trigger 1 and trigger 2, it means that
    % trigger 3 is 1 ms early. So we add +1
    if isempty(a) && isempty(c)
        trig2.ind(b) = trig2.ind(b) + 1;
        
    elseif isempty(a) && isempty(b)
        trig3.ind(c) = trig3.ind(c) + 1;
        
    elseif isempty(a) && ~isempty(b) && ~isempty(c)
        trig2.ind(b) = trig2.ind(b) + 1;
        trig3.ind(c) = trig3.ind(c) + 1;

    end
end

%% Re-compute triggers to check weird peaks again
triggerBinary = zeros(length(t), 3);
triggerBinary(trig1.ind,1) = 1;
triggerBinary(trig2.ind,2) = 1;
triggerBinary(trig3.ind,3) = 1;
% triggerBinary(trig3.ind,4) = 1;

trigger = triggerBinary * [1 2 4]'; clear triggerBinary;

%% Delete the blank triggers for now
for ii = 1:length(trigger)
    if trigger(ii) == 3 || trigger(ii) == 4;
        trigger(ii) = 0;
    end
end


%% check again at we have the right number of triggers
fprintf('%s\n','Check again for misaligned triggers')
fprintf('%s\t\t%s\n', 'Trigger nr', 'Sum')
for ii = 1:15; fprintf('%d\t%d\n', ii, sum(trigger == ii)); end


%% replot to make sure
if produce_figures

    figure(201); clf
    set(gca, 'ColorOrder', jet(4), 'Color', [.6 .6 .6]); hold all; 
    plot(bsxfun(@plus, ts(:,161:164), (1:4)*40000), '-o',  'LineWidth', 2); 
    legend(cellstr(num2str((161:164)')))

    weird = find(trigger == 9 | trigger == 10 | trigger == 11 | trigger == 12 | ...
    trigger == 13 | trigger == 14 | trigger == 15);
    set(gca, 'XTick', weird, 'XGrid', 'on')
    title('Trigger channels (indices and values)');
    xlabel('Time [s]');
    ylabel('Relative amplitude');
    legend('Trigger 1 / Chan 161', 'Trigger 2 / Chan 162', 'Trigger 3 / Chan 163', 'Trigger 4 / Chan 164');

end

%% Get epochs for every condition

% 432 triggers per condition = 
%   6 triggers per s * 6s / block * 2 blocks / experiment * 6 experiments

%% pad the time series because the last recording ended sligter before the last experiment ended

% Only for this run though
if strcmp(pwd,'/Volumes/server/Projects/MEG/SSMEG/02_SSMEG_02_28_2014');
    ts = padarray(ts, [2000, 0], NaN, 'post');
end

%% Get epoch indices from trigger values

if strcmp(experiment_name,'Attention')
    % Index of every trigger value
    trig_full_on_ind  = find(trigger == 1); 
    trig_left_on_ind  = find(trigger == 3);
    trig_right_on_ind = find(trigger == 5);   
elseif strcmp(experiment_name,'Default')
    % Index of every trigger value
    trig_full_on_ind  = find(trigger == 1); 
    trig_left_on_ind  = find(trigger == 5);
    trig_right_on_ind = find(trigger == 7);
    
end



cycles_per_s = 6; % six triggers in one epoch

% Get start value of every 1-s epoch
epoch_start_full_on_ind  = trig_full_on_ind(1:cycles_per_s:end); % So we have 72 epochs per condition
epoch_start_right_on_ind = trig_right_on_ind(1:cycles_per_s:end);
epoch_start_left_on_ind  = trig_left_on_ind(1:cycles_per_s:end);

ep = length(epoch_start_full_on_ind);

% Give each condition a number (Full = 1, Right = 2, Left = 3)
epoch_start_full_on_with_condition   = [epoch_start_full_on_ind, ones(length(epoch_start_full_on_ind),1), [1:1:ep]'];
epoch_start_right_on_with_condition  = [epoch_start_right_on_ind, 2*ones(length(epoch_start_right_on_ind),1), [ep+1:1:2*ep]'];
epoch_start_left_on_with_condition   = [epoch_start_left_on_ind, 3*ones(length(epoch_start_left_on_ind),1), [2*ep+1:1:3*ep]'];

% Define off conditions by on-condtions
num_epoch_time_pts = mode(diff(epoch_start_full_on_ind)); % frame duration is the same for every condition, i.e. 167 ms

epoch_start_full_off_ind  = epoch_start_full_on_ind + cycles_per_s * num_epoch_time_pts;
epoch_start_right_off_ind = epoch_start_right_on_ind + cycles_per_s * num_epoch_time_pts;
epoch_start_left_off_ind  = epoch_start_left_on_ind + cycles_per_s * num_epoch_time_pts;

% Give off conditions a number (Full off = 4, right off = 5, left off = 6)
epoch_start_full_off_with_condition   = [epoch_start_full_off_ind, 4*ones(length(epoch_start_full_off_ind),1), [3*ep+1:1:4*ep]'];
epoch_start_right_off_with_condition  = [epoch_start_right_off_ind, 5*ones(length(epoch_start_right_off_ind),1), [4*ep+1:1:5*ep]'];
epoch_start_left_off_with_condition   = [epoch_start_left_off_ind, 6*ones(length(epoch_start_left_off_ind),1), [5*ep+1:1:6*ep]'];

% % Concatenate the epochs with their condition numbers
epoch_list = vertcat(epoch_start_full_on_with_condition,epoch_start_right_on_with_condition, ...
             epoch_start_left_on_with_condition, epoch_start_full_off_with_condition, ...
                epoch_start_right_off_with_condition, epoch_start_left_off_with_condition);
         
% Sort them by time stamp
epochs_condition = sortrows(epoch_list,1);

% Save them
save epoch_conditions.mat epochs_condition

% Define for script
epochs_on = [epoch_start_full_on_ind, epoch_start_right_on_ind, epoch_start_left_on_ind];
epochs_off = [epoch_start_full_off_ind,epoch_start_right_off_ind,epoch_start_left_off_ind];



%% Check blank offsets
if produce_figures
    epoch_start_off_ind = sort(vertcat(epoch_start_full_off_ind,epoch_start_right_off_ind,epoch_start_left_off_ind));
    
    figure(202); clf
    set(gca, 'ColorOrder', jet(3), 'Color', [.6 .6 .6]); hold all; 
    plot(bsxfun(@plus, ts(:,161:163), (1:3)*40000), '-o',  'LineWidth', 2); 
    legend(cellstr(num2str((161:163)')))

    set(gca, 'XTick', epoch_start_off_ind, 'XGrid', 'on')
end

