%% s_GAMMA_MEG_visualize_natural
% visualize processed data from s_GAMMA_MEG_analysis.m
% TODO: loop across data sets

%% Parameters

% paths
if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end
project_pth    = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth       = '*_Gamma_*subj*';
% Find subject path
d = dir(fullfile(project_pth, data_pth));
%   restrict to directories
subj_pths = struct2cell(d);

% data parameters
fs                            = 1000;
intertrial_trigger_num        = 14;
which_session_to_visualize    = 17; %[7:12,14:16];
save_images                   = false;
using_denoised_data           = false;
suffix                        = 'preliminary';

for session_num = which_session_to_visualize
%% load data
condition_names = gamma_get_condition_names(session_num);
path_to_data = meg_gamma_get_path(session_num);
load_pth    = fullfile(path_to_data, 'processed');

if using_denoised_data
    save_pth = fullfile(project_pth,'/Images',subj_pths{session_num-1}, 'denoised');
    d        =  dir(fullfile(load_pth, '*_denoisedData_*boot*'));
    badChannels = [];
else
    save_pth = fullfile(project_pth,'/Images',subj_pths{session_num-1});
    d        =  dir(fullfile(load_pth, sprintf('*%s*', suffix)));
    badChannels = zeros(1,157);
end

results         = load(fullfile(load_pth, d(1).name)); % f x cond x chan

%% calculate contrasts 

num_conditions = size(results.w_pwr,2);
num_channels   = size(results.w_pwr,1);


% take the mean across channels 
if results.nboot > 1,
    summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
else
    summary_stat = @(x) nanmean(x,3);
end

end
