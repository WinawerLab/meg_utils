%% s_SSMEG_prepare_data
%
% This script loads sqd data from SSMEG experiment, epochs it, and resaves.
% The saved files from this script have the names:
%   s0*_conditions.mat
%   s0*_sensorData.mat
% Where s0* refers to each of subjects 1-8

%% Set analysis variables
project_pth                   = '/Volumes/server/Projects/MEG/SSMEG';
subjects                      = 1;
data_channels                 = 1:157;
trigger_channels              = 161:164;
photodiode_channel            = 192; 

fs                            = 1000;        % sample rate (Hz)
epoch_time                    = [0 1];       % start and end of epoch, relative to onset, in seconds
%% To run script, you need the Field trip toolbox

% Add fieldtrip path
field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
% meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
addpath(genpath(field_trip_pth))

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));

% Make sure there are exactly 8 subjects
assert(length(subject_pths)==8)

%% -------------------------------------
% ------- STRUCTURE THE RAW DATA -------
% --------------------------------------
%% Load in data
save_pth = fullfile(project_pth, 'manuscript_data');
if ~exist(save_pth, 'dir'), mkdir(save_pth); end

for which_subject = subjects
    data_pth = fullfile(project_pth, subject_pths(which_subject).name, 'raw');
    [ts, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG_*');
        
    if which_subject == 5        
        % This function is specifically made for session 8, look inside the
        % code if you want to use it for a different session!        
        onsets = ssmeg_get_triggers_from_photodiode(photodiode_channel, fs, ts);        
        
    else        
        trigger = meg_fix_triggers(ts(:,trigger_channels));  
        onsets = ssmeg_trigger_2_onsets(trigger, which_subject);
    end
    
    [sensorData, conditions] = meg_make_epochs(ts, onsets, epoch_time, fs);
    
    % Save behavior.mat file for source localization
    fullfield  = find(conditions==1);
    leftfield  = find(conditions==5);
    rightfield = find(conditions==7);
    blankfield = find(conditions==3);
    time.fullfield  = [0 1];
    time.leftfield  = [0 1];
    time.rightfield = [0 1];
    time.blankfield = [0 1];
    
%     save('/Volumes/server/Projects/MEG/SSMEG/04_SSMEG_04_01_2014_wl_subj006/behavior/behavior.mat', 'fullfield','leftfield','rightfield','blankfield','time');
    
%     save(fullfile(save_pth, sprintf('s%02d_sensorData', which_subject)), 'sensorData');
	
        
end