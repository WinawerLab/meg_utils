%% ssmeg_Analysis

% This script can be used to call MEG functions to analyse the data in several steps

%% Set analysis variables
project_pth                   = '/Volumes/server/Projects/MEG/SSMEG/fullOnly';

which_subjects                = 3;
denoise_with_nonphys_channels = false;       % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

produce_figures               = false;        % If you want figures in case of debugging, set to true
save_tseries                  = true;        % Save epoched time series?

denoise_via_pca               = false;        % Do you want to use PCA on a noise pool of channels to
experiment_name               = 'Default';   % Define condition of experiment, to define triggers (Choose between Attention or Default)

data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

fs                            = 1000;        % sample rate
epoch_time                    = [0 1];       % start and end of epoch, relative to trigger, in seconds

%% To run script, you need the Field trip toolbox

% Add fieldtrip path
% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip', 'yokogawa_defaults');

% Add SSMEG code path
addpath(fullfile(project_pth, 'code'));

% Find subjects for this project
subject_pths = dir(fullfile(project_pth, '*SSMEG_*'));

% Make sure there are exactly 8 subjects and one validation dataset
assert(length(subject_pths)==8+1)

%% -------------------------------------
% ------- STRUCTURE THE RAW DATA -------
% --------------------------------------
%% Load in data

data_pth = fullfile(project_pth, subject_pths(which_subjects).name, 'raw');

[ts, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG*');
%% Fix triggers

if which_subjects == 5
    pd_chan = 192; % Photodiode channel
    Fs = 1000; % Hz
    num_runs = 6; % Per single flicker/blank period
    num_epoch_time_pts = 1000;
    
    % This function is specifically made for session 8, look inside the
    % code if you want to use it for a different session!
    [trigger] = ssmeg_get_triggers_from_photodiode(pd_chan, Fs, num_runs, ts);
    
else
      
    trigger = meg_fix_triggers(ts(trigger_channels,:)');
    
end

onsets = ssmeg_trigger_2_onsets(trigger, 7, 'meg');
[sensorData, conditions] = meg_make_epochs(ts', onsets2, epoch_time, fs);

if save_tseries
    save(fullfile(project_pth, subject_pths(which_subjects).name, 'processed', sprintf('s%02d_sensorData.mat', which_subjects)), 'sensorData');
    save(fullfile(project_pth, subject_pths(which_subjects).name, 'processed', sprintf('s%02d_conditions', which_subjects)), 'conditions');    
end

return
%% -------------------------------------
% ------- PREPROCESS -------------------
% --------------------------------------

%% remove bad channels

%% Find bad epochs
if remove_bad_epochs
    
    % This identifies any epochs whos variance is outside some multiple of the
    % grand variance
    bad_epochs = meg_find_bad_epochs(sensorData(:,:,data_channels), [.05 20]);
    
    % any epoch in which more than 10% of channels were bad should be removed
    % entirely
    epochs_to_remove = mean(bad_epochs,2)>.1;
    
    % once we remove 'epochs_to_remove', check whether any channels have more
    % than 10% bad epochs, and we will remove these
    channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.1;
    
    bad_epochs(epochs_to_remove,:) = 1;
    bad_epochs(:,channels_to_remove) = 1;
    
    figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')
    
    sensorData = meg_remove_bad_epochs(bad_epochs, sensorData);
end

    
%% Denoise data by regressing out nuissance channel time series
if denoise_with_nonphys_channels      
      sensorData = meg_environmental_denoising(sensorData, environmental_channels,...
          data_channels, produce_figures);
end

%% DEAL WITH NANs, 1st epoch per block, etc

%% Save data for Helena's denoising
if save_tseries
    save(sprintf('s%02d_sensorData_preprocessed', which_subjects), 'sensorData');    
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






