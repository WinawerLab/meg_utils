% Get epochs from photo diode
function onsets = ssmeg_get_triggers_from_photodiode(pd_chan, Fs, ts)

cur_dir = pwd;
cd('/Volumes/server/Projects/MEG/SSMEG/08_SSMEG_06_20_2014_wl_subj011')

% Get length of run in time
t = (0:size(ts,1)-1) / Fs; % seconds

% Get peaks of photodiode response
data_pd   = diff(ts(:,pd_chan));
thresh_pd      = (max(data_pd) - min(data_pd)) / 50;
[pd.value_white, pd.ind_white] = findpeaks([0; -data_pd],'MINPEAKHEIGHT',thresh_pd);
[pd.value_black, pd.ind_black] = findpeaks([0; data_pd],'MINPEAKHEIGHT',thresh_pd);

pd.ind = sort([pd.ind_white; pd.ind_black]);

% Get rid of the runs where the photodiode triggers don't make sense
start_ind_1 = find(pd.ind == 488777);
pd.ind_shortened = pd.ind(start_ind_1:end);

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

conditions = [seq_3;seq_4;seq_5;seq_6;seq_7;seq_8;seq_9;seq_10;seq_11;seq_12;seq_13;seq_14;seq_15];

% % Give blank triggers a new nr, to match with the on conditions
% for ii = 1:length(master_seq)
%     if master_seq(ii) == 3 & master_seq(ii-1) == 1
%         master_seq(ii:ii+5) = 8; % full
%     elseif master_seq(ii) == 3 & master_seq(ii-1) == 5
%         master_seq(ii:ii+5) = 9; % right
%     elseif master_seq(ii) == 3 & master_seq(ii-1) == 7
%         master_seq(ii:ii+5) = 10; % left
%     end
% end

% Put triggers in same length timeseries
onsets = zeros(length(t), 1);
onsets(inds) = conditions; 

return

% Index of every trigger value
trig_full_on_ind  = find(trigger == 1); 
trig_left_on_ind  = find(trigger == 5);
trig_right_on_ind = find(trigger == 7);
trig_off_ind      = find(trigger == 3);

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

cd(cur_dir)
