%% BAIR MEG pilot analysis

% This script is a first attempt to analyse the MEG data in the BAIR
% project

% This script requires meg_utils and Fieldtrip


%% 0. Define general parameters

projectPth                    = '/Volumes/server/Projects/BAIR/MEG/data';
dataPth                       = 'wl_subj004_8.14.17';

dataChannels                  = 1:157;
environmentalChannels         = 158:160;
triggerChannels               = 161:168;
denoise_with_nonphys_channels = false;          % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;           % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true; 

fs                            = 1000;           % sample rate
epoch_start_end               = [0.00 1.00];  % start and end of epoch, relative to trigger, in seconds

%% 1. Import the data (raw sqd and behavioral files)

rawTs = meg_load_sqd_data(fullfile(projectPth, dataPth, 'raw'), '*Bair*');
rawTs = rawTs';

%% 2. Look at triggers

trigger = meg_fix_triggers(rawTs(:,triggerChannels));

%% 3. Make epochs

[ts, conditions]  = meg_make_epochs(rawTs, trigger, epoch_start_end, fs);



%% Preprocess (look at the triggers, remove bad channels, bad epochs)

% Plot some raw data
figure; for ii=1:157; plot(mean(ts(:,2603:2637,ii),2)); hold all; end
figure; for ii=1:157; plot(mean(ts(:,:,ii),2)); hold all; end
figure; for ii = 1:size(ts,2); plot(ts(:,ii,1)); hold all; end


tsClean = meg_environmental_denoising(ts);

figure; for ii=1:157; plot(mean(tsClean(:,2603:2637,ii),2)); hold all; end
figure; for ii=1:157; plot(mean(tsClean(:,:,ii),2)); hold all; end
figure; for ii = 1:size(ts,2); plot(tsClean(:,ii,1)); hold all; end


%% compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);

spectral_data = abs(fft(ts))/length(t)*2;
spectral_data_clean = abs(fft(tsClean))/length(t)*2;


figure; plot(f,mean(spectral_data(:,:,1),2)); xlim([1 200]); set(gca, 'XScale', 'log', 'YScale', 'log')
figure; plot(f,mean(spectral_data_clean(:,:,1),2)); xlim([1 200]); set(gca, 'XScale', 'log', 'YScale', 'log')

