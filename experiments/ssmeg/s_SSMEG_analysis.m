%% s_SSMEG_analysis

% This is a script to analyze SSMEG experiments, with use of the meg_utils
% functions


%%%%%%%%%%%%%%%%%%%%%%
%% Analysis options %%
%%%%%%%%%%%%%%%%%%%%%%
% Set analysis variables
project_pth                     = '/Volumes/server/Projects/MEG/SSMEG/';

% data to be analysed (% set nr)
data_pth                      = {'04_SSMEG_04_01_2014_wl_subj006', ... % 1
                                    '05_SSMEG_04_04_2014_wl_subj002', ... % 2
                                    '06_SSMEG_04_28_2014_wl_subj008', ... % 3
                                    '07_SSMEG_05_01_2014_wl_subj009', ... % 4
                                    '08_SSMEG_06_20_2014_wl_subj011', ... % 5
                                    '09_SSMEG_06_27_2014_wl_subj010', ... % 6 
                                    '10_SSMEG_08_12_2014_wl_subj004', ... % 7
                                    '11_SSMEG_08_13_2014_wl_subj005'}; % 8
num_data_sets                 = length(data_pth);

data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;
photodiode_channel            = 192;
    
denoise_with_nonphys_channels = false;       % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

produce_figures               = true;        % If you want figures in case of debugging, set to true

denoise_via_pca               = false;       % Do you want to use megdenoise?

fs                            = 1000;        % sample rate
% epoch_start_end               = [.05 0.55];  % start and end of epoch, relative to trigger, in seconds

save_images                   = false;

% condition names correspond to trigger numbers
condition_names               = {   'Full field' ... 
                                    'Left field' ...
                                    'Right field'};
                              
load_saved_epochs             = false;

which_data_sets_to_analyze = 1;


%% Add paths
if ~exist('meg_add_fieldtrip_paths.m')
    addpath(genpath('~/matlab/git/meg_utils'));
end

%change server-1 back to server
meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip', 'yokogawa_defaults')

%% Which subjects
for subject_num = which_data_sets_to_analyze
    
%%%%%%%%%%%%%%%%%%%
%% PREPROCESSING %%    
%%%%%%%%%%%%%%%%%%%
    
%% Make image folder
save_pth = fullfile(project_pth, 'Images', data_pth{subject_num});
if ~exist(save_pth, 'dir'), mkdir(session_image_pth); end


%% For new data sets
if load_saved_epochs == false;
    
    %% Load data (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, data_pth{subject_num}, 'raw'), '*SSMEG*');

    %% Extract triggers
    if subject_num == 5
        % Photodiode function --> generalize this function
    else
        trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    end

end

%% Instead of raw ts, one could also load in the preprocessed epochs



end