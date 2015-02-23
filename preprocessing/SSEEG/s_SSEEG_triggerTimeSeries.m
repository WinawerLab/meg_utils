
%% Script to extract DIN from NetStation event file and transform it into timeseries

% Note: Vistadisp code starts with a white square in the upper left corner
% to trigger photodiode.


%% Define parameters and load in event data
project_path = '/Volumes/server/Projects/EEG/SSEEG/';
session_name = 'Pilot_SSEEG_20150129_wl_subj001';
ev_data_file = 'Session_20150129_1007.evt';
el_data_file = 'Session_20150129_1007.mat';

pth = fullfile(project_path,'Data', session_name, 'raw', ev_data_file);
ev_time = eeg_read_evt_file(pth);

epoch_length     = 1;    % seconds
blocks_per_run   = 6;
trigs_per_block  = 36;
s_rate           = 1000; % Hz
first_run        = 2; % In case there are irrelevant runs recorderd to check stimulus code for presentation
trials_per_run   = 6;
blanks_per_trial = 6;                    
plot_figures     = true;

if plot_figures
    figure; plot(ev_time);
    xlabel('Nr of DIN')
    ylabel('Recorded time by NetStation')
end

%% Find startpoints of runs. Make trigger times increase continuously
% instead of resetting to zero at beginning of each run

startpoints = find((diff(ev_time) < 0))+1;
startpoints = [1 startpoints]; % no negative diff for the first run%% Make triggers increase continuously in time by accounting for zero reset
nr_runs = size(startpoints,2);

run_triggers = cell(nr_runs,1);
for ii = 1:nr_runs
    if ii == 1; % don't add anything
        run_triggers{ii} = ev_time(startpoints(ii):startpoints(ii+1)-1);
    elseif ii == length(startpoints) % add timepoints but take endpoint, since there is no next run
         run_triggers{ii} = sum(ev_time(startpoints(2:ii)-1)) + ev_time(startpoints(ii):end-1);
    else % add timepoints and go up to next run
         run_triggers{ii} = sum(ev_time(startpoints(2:ii)-1)) + ev_time(startpoints(ii):startpoints(ii+1)-1);
    end
end

these_inds = round(s_rate*cat(2,run_triggers{:}))+1;
clear run_triggers;

%% Load and concatenate electrode data
el_data = load(fullfile(project_path,'Data', session_name, 'raw', el_data_file));

cat_el_data = [];
fields = fieldnames(el_data);
for ii = 1:numel(fields)-3 % Assuming that there will be a DIN, TCPIP and Marks field
  cat_el_data{ii} = el_data.(fields{ii});
end

cat_el_data = cat(2,cat_el_data{:,:});
clear el_data;

%% create a time series with value of 100 for white, 0 for black square
num_time_pts = size(cat_el_data,2);
t = 0:num_time_pts-1;
ev_ts = zeros(1, num_time_pts);       
rf_rate = s_rate*mode(diff(ev_time)); % refresh rate in ms

for ii = 0:rf_rate
    ev_ts(these_inds+ii)=100; % 
end

if plot_figures
    figure; plot(t, ev_ts); hold on; plot(t,cat_el_data(1:4,:));
    xlabel('Time [ms]')
    ylabel('Electric field strength (microVolts and arbitrary units for triggers'); 
end

%% Cutoff irrelevant runs
if exist('first_run','var')
    preferred_start = ev_time(startpoints(first_run)-1)*1000;
    ev_ts = ev_ts(preferred_start:end);
    t = t(preferred_start:end);
    cat_el_data_short = cat_el_data(:,preferred_start:end);
end 

if plot_figures 
        figure; plot(t, ev_ts); hold on; plot(t,cat_el_data_short(1:4,:));
            xlabel('Time [ms]')
            ylabel('Electric field strength (microVolts and arbitrary units for triggers'); 
end 

%% Locate initializing sequence
ev_inds = find(ev_ts == 100); 

% Find all indices that have a difference 1 ms higher than the rf_rate
higher_than_rf = ev_inds(diff(ev_inds)>(rf_rate+1)); 

% We will take all the indices for whom rf_rate < diff < 80 ms
unique_diff = unique(diff(t(higher_than_rf)));
diff_DIN_80 = unique_diff<80;

% For this data set the first 4 indices, time differences are 49,50,66,67 
% and will be used to get the time of first DIN's of each run's initiation 
% sequence.
init_seq_diff_ind = [];
for ii = 1:sum(diff_DIN_80);
    init_seq_diff_ind{ii} = find(diff(t(higher_than_rf))==unique_diff(ii));
end

init_seq_diff_ind = sort(cell2mat(init_seq_diff_ind));

% Since we only need the first, and the sequence will always be 4 flickers
end_init_seq = init_seq_diff_ind(3:3:end);

%% Define the start of each run as 
start_run_ind  = higher_than_rf(end_init_seq+9);
start_run_time = t(start_run_ind);

triggers_start_time = [];
triggers_start_ind = [];
for ii = 1:size(start_run_time,2)
    temp = find(t(ev_inds) == start_run_time(ii));
    
    triggers_start_ind{ii} = ev_inds(temp+1);
    triggers_start_time{ii} = t(ev_inds(temp+1));
    
end

triggers_start_ind = cell2mat(triggers_start_ind);
triggers_start_time = cell2mat(triggers_start_time);
clear start_run_ind    start_run_time    end_init_seq    init_seq_diff_ind; 

%% 


%% potentially useful still
% 
% t0    =  (sum(ev_time(startpoints(2:first_run)-1))*1000) + 1;
% t_end =  size(cat_el_data,2);
% a = ev_ts(t0)
% ev_ts =  zeros(1,(t_end - t0));

% Cutoff irrelevant runs
% if exist('first_run','var')
%     preferred_start = ev_time(startpoints(first_run)-1)*1000;
%     ev_ts = ev_ts(preferred_start:end);
%     t = t(preferred_start:end);
%     
%     if plot_figures
%         % Plot only relevant data
%         figure; plot(t, ev_ts); hold on; plot(t,cat_el_data(1:4,preferred_start:end));
%         xlabel('Time [ms]')
%         ylabel('Electric field strength (microVolts and arbitrary units for triggers'); 
%     end
%     
% end

%% create new time series which begins at t=0
% t2 = linspace(1,size(t,2),size(t,2));       
%     if plot_figures
%         % Plot ev_ts and a few visual channels along t2
%         figure; plot(t2, ev_ts); hold on; plot(t2,cat_el_data_short(69:71,:));
%         hold on; plot(triggers_start_ind,100,'ro');
%             xlabel('Time [ms]')
%             ylabel('Electric field strength (microVolts and arbitrary units for triggers'); 
%     end
% 
% ev_inds_start = zeros(1,size(triggers_start_ind,2));
% for ii = 1:size(triggers_start_ind,2)
%         ev_inds_start(1,ii) = find(ev_inds == triggers_start_ind(ii));
% end

%% define triggers as every fifth DIN, and define Epochs as every sixth trigger
% 
% num_runs = size(ev_inds_start,2);
% cycles_per_run = (trigs_per_block * blocks_per_run);
% triggers_all = zeros(triggers_per_run,num_runs);
% blanks_per_run = (trials_per_run * blanks_per_trial);
% blank_start_times = zeros(blanks_per_run,num_runs);
% cycle_length = 1/Hz;
% DINS_per_cycle = round(cycle_length/rf_rate/2);
% 
% for ll = 1:num_runs
%         jj = 1;
%             for ii = ev_inds_start(ll):ev_inds_start(ll)+(cycles_per_run*DINS_per_cycle)-1
%                 if (ev_inds(ii) - ev_inds(ii-1)) > (rf_rate*1000+1) 
%                 triggers_all(jj,ll) = ev_inds(ii);
%                     jj = jj+1;
%                 end
%             end        
% end            
% 
% if plot_figures
%     for ii = 1:num_runs
%     plot(triggers_all(:,ii), 99, 'rx');
%     end
% end
% 
% epochs = zeros(cycles_per_run/6,num_runs);
% for ii = 1:num_runs
%     epochs(:,ii) = triggers_all(1:6:end,ii);
% end
% 
% if plot_figures
%     for ii = 1:num_runs
%     plot(epochs(:,ii), 98, 'kx');
%     end
% end

%% create a 
% epoched_data = zeros(numel(epochs)-2,epoch_length*1000);
% 
% for ll = 1:size(epoched_data,1)
%         for ii = epoched_data(
%             epoched_data(ll,ii) = cat_el_data_short(74,epochs();
% end
%     
%% probably trash

%(ev_inds(ii) - ev_inds(ii-1)) < epoch_length*2*1000
%for ll = 1:num_runs
%    if ll < num_runs
%         jj = 1;
%             for ii = ev_inds_start(ll):ev_inds_start(ll)+triggers_per_run
%                 if (ev_inds(ii) - ev_inds(ii-1)) > (rf_rate*1000+1) 
%                 triggers_all(jj,ll) = ev_inds(ii);
%                     jj = jj+1;
%                 end
%             end
%     elseif ll == num_runs 
%         jj = 1;
%             for ii = ev_inds_start(num_runs):
%                 if (ev_inds(ii) - ev_inds(ii-1)) > (rf_rate*1000+1) 
%                 triggers_all(jj,ll) = ev_inds(ii);
%                     jj = jj+1;
%                 end
%             end
%     end
% end

% plot(triggers_start_time,zeros(length(triggers_start_time),1)`,'ro');
% plot(t([ev_inds(2064):16:ev_inds(2064+(16*36))]),zeros(length(t([ev_inds(2064):5:ev_inds(2064+(16*36))])),1),'go')

% all_trigger_inds = ev_ts(ev_inds(triggers_start_ind(1)));
% :4:triggers_start_ind(1)+36)
% triggers = triggers_start_epoch:5:triggers_start_epoch+36)