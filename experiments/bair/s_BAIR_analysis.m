%% BAIR MEG pilot analysis

% This script is a first attempt to analyse the MEG data in the BAIR
% project


% NB. This script requires meg_utils and Fieldtrip


% Info about the experiment:

% HRF (zebra)
%  Stimuli are presented for 0.25 seconds, with an exponentially
%  distributed interstimulus interval (mean ~9s, range 3-24s).
% HRF
%       Exponential ISI, mean ~9s, range 3-24s
% TASK
%       12 s fixation task, 12 s relax, repeat 8 times (to measure anticipatory BOLD)
% RETINOTOPY
%       Like Dumoulin and Wandell, 2008: 8 sweeps ? 4 cardinal, 4 diagonal (diagonals include 50% blanks)
% SPATIOTEMPORAL
%        VISUAL: 36 unique stimuli, shown once each per scan (0.5 s except for temporal stimuli),
%                with mean ISI of 4.5 s, range 3-6 s; orientation (3; 1 grating, 1 plaid, 1 circular);
%                contrast (5; noise patterns); spacing: (5: noise patterns, 1 overlaps with contrast);
%                objects (12: 3 faces, 3 scenes, 3 objects, 3 bodies); temporal (12; 6 durations; 6 ISIs);

% SPATIO TEMPORAL (12 repeats, ECoG, fMRI; 24 for E/MEG !?)
% For visual experiments, we use band-pass, gray-scale images, spanning
% many stimulus dimensions. Twelve were used in a prior publication [69,70],
% varying in contrast, number of component orientations (1, 2 or 16
% superimposed gratings), or spacing between contrast elements (from very
% sparse to very dense). Twelve are natural images of faces, objects, and
% scenes (also gray-scale, band-pass). These stimuli will be presented for
% 0.5 seconds each. Twelve other stimuli are simple noise patterns shown
% with different temporal profiles (single pulses with variable duration;
% or multiple pulses with variable interstimulus interval).


% CRF      - 5 (zebra)                     KNK 162 164 166 168 116
% Orient   - 3 (grating, plaid, circular)  KNK 150, 154, 158 (*HC)
% Sparsity - 4 (zebras)                    KNK 181 182 183 184
% 1 Pulse  - 6 (zebra??)                   KNK 183 * 6
% 2 Pulses - 6 (zebra??)                   KNK 183 * 6

% Faces -    4                             KNK 171 (sample 6 * 8 for 12 runs, 4 each)
% Letters -  4                             KNK 173 (sample 6 * 8 for 12 runs, 4 each)
% Scenes -   4                             KNK 175 (sample 6 * 8 for 12 runs, 4 each)

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
rawTs = rawTs'; % transpose data because meg_fix_trigges requires time x channels

%% 2. Look at triggers

trigger = meg_fix_triggers(rawTs(:,triggerChannels));

%% 3. Make epochs

% [UNDER CONSTRUCTION]: Think about how to make epochs (deal with different
% experiments and epoch lengths)
[ts, conditions]  = meg_make_epochs(rawTs, trigger, epoch_start_end, fs);


%% Preprocess (look at the triggers, remove bad channels, bad epochs)

% Plot some raw data
figure; for ii=1:157; plot(mean(ts(:,2603:2637,ii),2)); hold all; end 
figure; for ii=1:157; plot(mean(ts(:,:,ii),2)); hold all; end
figure; for ii = 1:size(ts,2); plot(ts(:,ii,1)); hold all; end

% Clean up data by projecting out environmental channels
tsClean = meg_environmental_denoising(ts);

% plot some clean data for comparison
figure; for ii=1:157; plot(mean(tsClean(:,2603:2637,ii),2)); hold all; end
figure; for ii=1:157; plot(mean(tsClean(:,:,ii),2)); hold all; end
figure; for ii = 1:size(ts,2); plot(tsClean(:,ii,1)); hold all; end


%% compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);

spectral_data = abs(fft(ts))/length(t)*2;
spectral_data_clean = abs(fft(tsClean))/length(t)*2;

figure; plot(f,mean(spectral_data(:,:,1),2)); xlim([1 200]); set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Frequency (Hz)'); ylabel('Field strength (Tesla)')
figure; plot(f,mean(spectral_data_clean(:,:,1),2)); xlim([1 200]); set(gca, 'XScale', 'log', 'YScale', 'log'); 
xlabel('Frequency (Hz)'); ylabel('Field strength (Tesla)')

%% Behavioral files

% load one of the behavioral files.
d = dir(fullfile(projectPth, dataPth, 'behavior'));

behavior_file = load(fullfile(projectPth, dataPth, 'behavior', d(3).name));