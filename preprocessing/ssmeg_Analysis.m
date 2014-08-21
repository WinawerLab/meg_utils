%% ssmeg_Analysis

% This script can be used to call MEG functions to analyse the data in several steps

%% Set analysis variables

denoise_with_nonphys_channels = false;       % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

produce_figures               = true;        % If you want figures in case of debugging, set to true
save_tseries                  = true;        % Save epoched time series?

denoise_via_pca               = false;        % Do you want to use PCA on a noise pool of channels to
% denoise the signal
experiment_name               = 'Default';   % Define condition of experiment, to define triggers (Choose between Attention or Default)

%% To run script, you need the Field trip toolbox

% Add fieldtrip path
if isempty(which('ft_analysispipeline')),
    addpath('/Volumes/server/Projects/MEG/code/fieldtrip');
    
    % Add the other fieldtrip paths that are not automatically added:
    addpath('/Volumes/server/Projects/MEG/code/fieldtrip/external/yokogawa');
    addpath('/Volumes/server/Projects/MEG/code/fieldtrip/external/sqdproject');
    addpath('/Volumes/server/Projects/MEG/code/fieldtrip/fieldtrip_private/');
    addpath('/Volumes/server/Projects/MEG/code/fieldtrip/compat/');
    
    % and then run the ft_defaults which will load the necessary paths:
    ft_defaults
end

% Add SSMEG code path
addpath('/Volumes/server/Projects/MEG/SSMEG/code');

%% **** Preprocessing ****

%% Load in data
data_pth = '/Volumes/server/Projects/MEG/SSMEG/11_SSMEG_08_13_2014_wl_subj005/';
[ts, meg_files, stim_file, num_runs] = ssmeg_loaddata(data_pth);

%% Fix triggers

if strcmp(data_pth, '/Volumes/server/Projects/MEG/SSMEG/08_SSMEG_06_20_2014_subj011/')
    pd_chan = 192; % Photodiode channel
    Fs = 1000; % Hz
    num_runs = 6; % Per single flicker/blank period
    num_epoch_time_pts = 1000;
    
    % This function is specifically made for session 8, look inside the
    % code if you want to use it for a different session!
    [t, epochs_on,epochs_off] = ssmeg_get_triggers_from_photodiode(pd_chan, Fs, num_runs, ts);
    
else
%     [t, ts, num_epoch_time_pts,epochs_on,epochs_off,epochs_condition] = ssmeg_fix_triggers(...
%         ts, experiment_name, produce_figures);
      
    [t, ts, num_epoch_time_pts,epochs_on,epochs_off,epochs_condition] = ssmeg_fix_triggers_v2(...
        ts, experiment_name, produce_figures);

end

%% Make epochs

[ts_on, ts_off, num_epochs] = ssmeg_make_epochs(ts, num_epoch_time_pts,epochs_on,epochs_off);
%% remove bad channels

if remove_bad_channels
    bad_channels = ssmeg_find_bad_channels(ts_on(:,:,1:157));
    ts_on  = ssmeg_remove_bad_channels(bad_channels, ts_on);
    ts_off = ssmeg_remove_bad_channels(bad_channels, ts_off);
end

%% Denoise data by regressing out nuissance channel time series

% Denoise data with 3 noise channels
if denoise_with_nonphys_channels
    if ~exist('denoised_with_nuissance_data.mat', 'file')
        [ts_on, ts_off] = ssmeg_denoising(ts_on,ts_off, num_epochs, produce_figures);
    else fprintf('Loading data.. This may take a couple of seconds\n');
        load(fullfile(data_pth,'denoised_with_nuissance_data.mat'));
    end
end

%% Remove epochs with outliers

% Make matrices with marking the 'bad' epochs.
% Use variance to remove epochs with outliers.

if remove_bad_epochs
    [bad_epochs_on,bad_epochs_off,ts_on,ts_off] = ...
        ssmeg_bad_epochs(t, num_epoch_time_pts, ts_on, ts_off, num_epochs, produce_figures);
end


%% Check for bad channels again:

if remove_bad_channels
    bad_channels = ssmeg_find_bad_channels(ts_on(:,:,1:157));
    bad_channels = [bad_channels];
    ts_on  = ssmeg_remove_bad_channels(bad_channels, ts_on);
    ts_off = ssmeg_remove_bad_channels(bad_channels, ts_off);
end

%% Save data for Helena's denoising
if save_tseries
    
    ts_on_full   = ts_on(:,1:num_epochs,:); % 1
    ts_on_right  = ts_on(:,(num_epochs+1):2*num_epochs,:); % 2
    ts_on_left   = ts_on(:,(2*num_epochs)+1:3*num_epochs,:); % 3
    
    ts_off_full   = ts_off(:,1:num_epochs,:); % 4
    ts_off_right  = ts_off(:,(num_epochs+1):2*num_epochs,:); % 5
    ts_off_left   = ts_off(:,(2*num_epochs)+1:3*num_epochs,:); % 6
    
     save data_no_nonphys_denoising.mat ts_on_full ...
         ts_on_right ts_on_left ...
         ts_off_full ts_off_right ts_off_left
    
    
    save ts_on_full.mat ts_on_full
    save ts_on_left.mat ts_on_left
    save ts_on_right.mat ts_on_right
    
    save ts_off_full.mat ts_off_full
    save ts_off_left.mat ts_off_left
    save ts_off_right.mat ts_off_right
    
    
end

%% In case you want to lose first and last epoch of every flicker period
keep_epochs = true(1,num_epochs);
keep_epochs(1:6:end) = false;
keep_epochs(6:6:end) = false;

%% Fourier analysis

okEpochs = [];

[amps_on_full,amps_on_right,amps_on_left, ...
    amps_off_full,amps_off_right,amps_off_left] = ...
    ssmeg_fourier(t, num_epoch_time_pts, ts_on, ts_off, num_epochs, keep_epochs, okEpochs, produce_figures);

%% Calculate coherence

plotType = '2d'; % any of 2d, 3d, both
produceFigures = true;
[coh_full, coh_right, coh_left] = ssmeg_cal_coh(...
    amps_on_full,amps_on_right,amps_on_left, meg_files, plotType, keep_epochs, produceFigures);

%% meg_denoise?

if denoise_via_pca
    
%     sl_freq = 12; % Hz
%     design_nrs_all = [1 2 3 4 5 6]; % Use three on and three off conditions
%                                     % In case you want one condition, for
%                                     % example full field, use [1 4];
% 
%     % Do the MEG Denoise
%     [results,evalout,denoisedspec,denoisedts, badChannels, okEpochs] = meg_denoise_pca(num_epoch_time_pts, sl_freq, design_nrs_all);
%     
%     % Reshape the matrix so we get the old timepoints x epochs x channel
%     % structure back.
%     denoised_ts_all = permute(denoisedts{1},[2,3,1]);
%     
%     % Look which epochs are removed in the denoising, so we can redefine
%     % conditions
%     bad_epochs = find(okEpochs == 0);
%     
%     nconditions = 6;
%     nr_pre_epochs = num_epochs*nconditions;
%     halfway = nr_pre_epochs/2;
%     
%     % First define on an off conditions
%     bad_on_epochs = 0;
%     bad_off_epochs = 0;
%     for ii = 1:length(bad_epochs)
%         if bad_epochs(ii) < halfway+1
%             bad_on_epochs = bad_on_epochs +1;
%         else
%             bad_off_epochs = bad_off_epochs +1;
%         end
%     end
%     
%     ts_denoised_on = denoised_ts_all(:,1:halfway-bad_on_epochs,:);
%     ts_denoised_off = denoised_ts_all(:,(halfway-bad_on_epochs):end,:);
%     
%     % Find removed_channels
%     removed_channels = find(badChannels == 1);
%     
%     % Predefine whole timeseries matrix
%     ts_on = zeros(num_epoch_time_pts,size(ts_denoised_on,2),157);
%     ts_off = zeros(num_epoch_time_pts,size(ts_denoised_off,2),157);
%     
%     % Add NaN's for removed channels
%     chan_nr = 1;
%     for ii = 1:length(badChannels)
%         if badChannels(ii) == 0
%             ts_on(:,:,ii) = ts_denoised_on(:,:,chan_nr);
%             ts_off(:,:,ii) = ts_denoised_off(:,:,chan_nr);
%             chan_nr =  chan_nr +1;
%             
%         elseif badChannels(ii) == 1;
%             ts_on(:,:,ii) = NaN;
%             ts_off(:,:,ii) = NaN;
%         end
%     end
%     
%     
%     % Re-calculate the amplitudes in the spectrum after denosing
%     [amps_on_full,amps_on_right,amps_on_left, ...
%         amps_off_full,amps_off_right,amps_off_left] = ...
%         ssmeg_fourier(t, num_epoch_time_pts, ts_on, ts_off, num_epochs, keep_epochs, okEpochs, produce_figures);
end

%% Redo coherence for every principle component

plotType = '2d'; % any of 2d, 3d, both
produceFigures = true;
[coh_full, coh_right, coh_left] = ssmeg_cal_coh(...
    amps_on_full,amps_on_right,amps_on_left, meg_files, plotType, keep_epochs, produceFigures);

%% Compute broadband & spectra
% amps = cat(4, amps_on_full(:,keep_epochs(:),:), amps_off_full(:,keep_epochs(:),:));

% Check the epoch size for every condition and add NaN's, so we can use the
% same function for denoised and non-denoised data.
[amps_on_full, amps_off_full, ...
    amps_on_right, amps_off_right, ...
    amps_on_left, amps_off_left] = ssmeg_check_epoch_sizes(num_epochs,amps_on_full, amps_off_full, ...
    amps_on_right, amps_off_right, ...
    amps_on_left, amps_off_left);

amps = cat(4, amps_on_full, amps_off_full);
[ab, log_frequency] = ssmeg_broadband(num_epoch_time_pts, amps);
ssm_plotOnMesh(-diff(ab)./(sum(ab)),'full field broadband diff/sum', 300, [], plotType);
set(gca,'CLim',[-0.1 0.1]);

amps = cat(4, amps_on_left(:,:,:), amps_off_left(:,:,:));
[ab, log_frequency] = ssmeg_broadband(num_epoch_time_pts, amps);
ssm_plotOnMesh(-diff(ab)./(sum(ab)),'left frield broadband diff/sum', 301, [], plotType);
set(gca,'CLim',[-0.1 0.1]);

amps = cat(4, amps_on_right(:,:,:), amps_off_right(:,:,:));
[ab, log_frequency] = ssmeg_broadband(num_epoch_time_pts, amps);
ssm_plotOnMesh(-diff(ab)./(sum(ab)),'right field broadband diff/sum', 302, [], plotType);
set(gca,'CLim',[-0.1 0.1]);

%% channel loop
amps = cat(4, amps_on_full(:,keep_epochs(:),:), amps_off_full(:,keep_epochs(:),:));
amps = cat(4, amps_on_full, amps_off_full);

[ab, log_frequency, bbPoly, frequencies] = ssmeg_broadband(num_epoch_time_pts, amps);

sensor_data = -diff(ab)./(sum(ab));
npcs=0;
ssmeg_make_spectra_figures(sensor_data, data_pth, ab, bbPoly, frequencies, amps_on_full, amps_off_full, npcs, denoise_via_pca)






