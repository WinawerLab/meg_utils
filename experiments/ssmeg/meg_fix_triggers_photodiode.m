function trigger = meg_fix_triggers_photodiode(raw_ts, start_ind, threshold, fs, photodiode_channel)

%% DON'T USE FOR NOW, ITS NOT READY YET

%% Description of meg_fix_triggers_photodiode

% trigger = meg_fix_triggers_photodiode(raw_ts, start_ind, [threshold], [fs], [photodiode_channel])

% This function is based on the function meg_fix_triggers, but now you define 
% triggers by the use of photodiode timeseries. For now, one needs to input 
% the raw_data and a start index. That is, because the continuous recording 
% starts before the experiments was started, there may be some unrelated 
% 'triggers' in the timeseries that may vary per experiment.

% INPUTS:
% raw_ts    =   Timeseries [timepoints x channels]
% start_ind =   one integer that defines your starting point of the first
%               run.

% OUTPUTS: 





%% Deal with inputs
if ~exist(threshold, 'var') || isempty(threshold), threshold = 50; end % threshold to define peak 
if ~exist(fs, 'var') || isempty(fs), fs = 1000; end % sampling rate
if ~exist(photodiode_channel, 'var') || isempty(photodiode_channel), photodiode_channel = 192; end % sampling rate


%% Get length of run in time

%% OLD CODE
% t = (0:size(raw_ts,1)-1) / fs; % seconds

% % Get peaks of photodiode response
% data_pd             = diff(raw_ts(:,photodiode_channel));
% thresh_data_pd      = (max(data_pd) - min(data_pd)) / threshold;
% [pd.value_white, pd.ind_white] = findpeaks([0; -data_pd],'MINPEAKHEIGHT',thresh_data_pd);
% [pd.value_black, pd.ind_black] = findpeaks([0; data_pd],'MINPEAKHEIGHT',thresh_data_pd);
% 
% % Sort the triggers in time
% pd.ind = sort([pd.ind_white; pd.ind_black]);
% 
% % Get rid of the runs where the photodiode triggers don't make sense
% start_value = pd.ind(start_ind);
% pd.ind_shortened = pd.ind(start_ind:end);
% 
% % differentiate to isolate trigger onsets, and pad to preserve length
% ts = padarray(diff(pd.ind_shortened), [1 0], 0, 'post');

%%

ts = raw_ts(:,photodiode_channel);

% rescale to [0 1]
ts = ts - min(ts(:));
ts = ts / max(ts(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(ts)) == 0, trigger_is_high = true; 
else                     trigger_is_high = false; end 

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, ts = 1 - ts; end

% threshold to binarize from analog recording 
ts = ts > 0.1;

% differentiate to isolate trigger onsets, and pad to preserve length
ts = padarray(diff(ts), [1 0], 0, 'post');

% rectify 
ts(ts<0) = 0;

% mark the samples as 1 when any trigger channel signaled 
any_trigger      = sum(ts,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);

% Delete the stuff you don't want *This needs to be implemented in the code
% some day*
any_trigger_inds = any_trigger_inds(start_ind:end);

% check for really short triggers
assert(sum(diff(any_trigger_inds) < 10)==0)

%%

% get triggers around mode of trigger indices
criterion =  mode(diff(any_trigger_inds)) + [-2:2];
triggers_to_keep = any_trigger_inds==criterion;

difference_inds = diff(any_trigger_inds);
% for ii = 1:12
%     if ii == 1
%         inds.on(ii,:) = pd.ind_shortened(1:12:(12*36)-1);
%     else
%         inds.on(ii,:) = pd.ind_shortened(1:12:(12*36)-1);
%     end

% Define on periods (To Do: Automate this)
run3_on = pd.ind_shortened(1:12:431);

run4_on = pd.ind_shortened(441:12:872);

run5_on = pd.ind_shortened(881:12:1312);

run6_on = pd.ind_shortened(1321:12:1742);

run7_on = pd.ind_shortened(1761:12:2184);

run8_on = pd.ind_shortened(2201:12:2622);

run9_on = pd.ind_shortened(2641:12:3062);

run10_on = pd.ind_shortened(3081:12:3503);

run11_on = pd.ind_shortened(3521:12:3942);

run12_on = pd.ind_shortened(3961:12:4389);

run13_on = pd.ind_shortened(4401:12:4822);

run14_on = pd.ind_shortened(4841:12:5262);

run15_on = pd.ind_shortened(5281:12:5702);

% Concatenate indices and define off periods
nr_cycles = 6;
nr_epoch_timepoints = 1000;

inds_on = [run3_on; run4_on; run5_on; run6_on; run7_on; run8_on; run9_on; run10_on; run11_on; run12_on; run13_on; run14_on; run15_on];
inds_off = inds_on + nr_cycles*nr_epoch_timepoints;

inds = sort([inds_on;inds_off]);

% Plot all triggers 
figure(200); clf
set(gca, 'ColorOrder', jet(3), 'Color', [.6 .6 .6]); hold all;
plot(bsxfun(@plus, ts(:,192), (1)*40000), '-o',  'LineWidth', 2);
legend(cellstr(num2str((192)')))

set(gca, 'XTick', inds, 'XGrid', 'on');
title('Trigger channels (indices and values)');
xlabel('Time [s]');
ylabel('Relative amplitude');

% Get the stimulus sequence from the behavior mat-files.
mat_files = dir('behavior/*.mat');
stim_file = [];
stim_file.run3 = load(fullfile('behavior',mat_files(3).name)); 
stim_file.run4 = load(fullfile('behavior',mat_files(4).name)); 
stim_file.run5 = load(fullfile('behavior',mat_files(5).name)); 
stim_file.run6 = load(fullfile('behavior',mat_files(6).name)); 
stim_file.run7 = load(fullfile('behavior',mat_files(7).name)); 
stim_file.run8 = load(fullfile('behavior',mat_files(8).name)); 
stim_file.run9 = load(fullfile('behavior',mat_files(9).name)); 
stim_file.run10 = load(fullfile('behavior',mat_files(10).name)); 
stim_file.run11 = load(fullfile('behavior',mat_files(11).name)); 
stim_file.run12 = load(fullfile('behavior',mat_files(12).name)); 
stim_file.run13 = load(fullfile('behavior',mat_files(13).name));
stim_file.run14 = load(fullfile('behavior',mat_files(14).name));
stim_file.run15 = load(fullfile('behavior',mat_files(15).name));

% Load sequences. % 1,2 = full / 3,4 = blank / 5,6 = Left / 7,8 = right
seq_3 = stim_file.run3.stimulus.seq(1:12:end);
seq_4 = stim_file.run4.stimulus.seq(1:12:end);
seq_5 = stim_file.run5.stimulus.seq(1:12:end);
seq_6 = stim_file.run6.stimulus.seq(1:12:end);
seq_7 = stim_file.run7.stimulus.seq(1:12:end);
seq_8 = stim_file.run8.stimulus.seq(1:12:end);
seq_9 = stim_file.run9.stimulus.seq(1:12:end);
seq_10 = stim_file.run10.stimulus.seq(1:12:end);
seq_11 = stim_file.run11.stimulus.seq(1:12:end);
seq_12 = stim_file.run12.stimulus.seq(1:12:end);
seq_13 = stim_file.run13.stimulus.seq(1:12:end);
seq_14 = stim_file.run14.stimulus.seq(1:12:end);
seq_15 = stim_file.run15.stimulus.seq(1:12:end);

master_seq = [seq_3;seq_4;seq_5;seq_6;seq_7;seq_8;seq_9;seq_10;seq_11;seq_12;seq_13;seq_14;seq_15];

% Give blank triggers a new nr, to match with the on conditions
for ii = 1:length(master_seq)
    if master_seq(ii) == 3 & master_seq(ii-1) == 1
        master_seq(ii:ii+5) = 8; % full
    elseif master_seq(ii) == 3 & master_seq(ii-1) == 5
        master_seq(ii:ii+5) = 9; % right
    elseif master_seq(ii) == 3 & master_seq(ii-1) == 7
        master_seq(ii:ii+5) = 10; % left
    end
end

% Put triggers in same length timeseries
trigger = zeros(length(t), 1);
trigger(inds) = master_seq; 

% Index of every trigger value
trig_full_on_ind  = find(trigger == 1); 
trig_left_on_ind  = find(trigger == 5);
trig_right_on_ind = find(trigger == 7);

% Rename and define
epoch_start_full_on_ind  = trig_full_on_ind; 
epoch_start_right_on_ind = trig_right_on_ind;
epoch_start_left_on_ind  = trig_left_on_ind;
        
epoch_start_full_off_ind = find(trigger == 8);
epoch_start_right_off_ind = find(trigger == 9);
epoch_start_left_off_ind = find(trigger == 10);

ep = length(epoch_start_full_on_ind);

% Give each condition a number (Full = 1, Right = 2, Left = 3)
epoch_start_full_on_with_condition   = [epoch_start_full_on_ind, ones(length(epoch_start_full_on_ind),1), [1:1:ep]'];
epoch_start_right_on_with_condition  = [epoch_start_right_on_ind, 2*ones(length(epoch_start_right_on_ind),1), [ep+1:1:2*ep]'];
epoch_start_left_on_with_condition   = [epoch_start_left_on_ind, 3*ones(length(epoch_start_left_on_ind),1), [2*ep+1:1:3*ep]'];

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

end


