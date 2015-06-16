% s_GAMMA_MEG_evoked_responses

% Analyze and visualize evoked responses from MEG Gamma experiments. Subjects saw
% several kinds of stimuli, including gratings of various spatial
% frequencies, plaids, and noise patterns (phase-randomized with 1/f^n
% spectral power distrubutions)
%
% Stimuli were static, and on the screen either for 500 ms (subjects 1-3)
% or 1000 ms (subjects 4-6) with 500 ms ISI.
%

% Analysis options
%% Set analysis variables
project_pth                     = '/Volumes/server/Projects/MEG/Gamma/Data';

% data to be analysed
data_pth                      = '*_Gamma_*subj*';

data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

denoise_with_nonphys_channels = true;         % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;         % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;         % Remove channels whose median sd is outside some range

nboot                         = 1;            % number of bootstrap samples

produce_figures               = true;         % If you want figures in case of debugging, set to true

fs                            = 1000;         % sample rate
epoch_start_end               = [0.000 0.999];% start and end of epoch, relative to trigger, in seconds

intertrial_trigger_num        = 11;           % the MEG trigger value that corresponds to the intertrial interval

save_images                   = false;

% condition names correspond to trigger numbers
condition_names               = {   ...
    'White Noise' ...
    'Binarized White Noise' ...
    'Pink Noise' ...
    'Brown Noise' ...
    'Gratings(0.36 cpd)' ...
    'Gratings(0.73 cpd)' ...
    'Gratings(1.45 cpd)' ...
    'Gratings(2.90 cpd)' ...
    'Plaid(1.45 cpd)'...
    'Blank'};

which_data_sets_to_analyze = 5;
blank_condition = strcmpi(condition_names, 'blank');

%% Add paths

%change server-1 back to server
if isempty(which('ft_defaults'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',...
        {'yokogawa', 'sqdproject'})
end

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

%% Loops over datasets
for subject_num = which_data_sets_to_analyze
    
    save_pth = fullfile(project_pth, 'Images', subj_pths{subject_num});
    if ~exist(save_pth, 'dir'), mkdir(save_pth); end
    
    % --------------------------------------------------------------------
    % ------------------ PREPROCESS THE DATA -----------------------------
    % --------------------------------------------------------------------
    %% Load data (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject_num}, 'raw'), '*Gamma*');
    
    %% Extract triggers
    trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    
    %% Make epochs
    [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);        
    % remove intertrial intervals
    iti               = conditions == intertrial_trigger_num;
    ts                = ts(:,~iti, :);
    conditions        = conditions(~iti);
    conditions_unique = unique(conditions);
    num_conditions    = length(condition_names);
    num_time_points   = size(ts, 1);
    %% Remove bad epochs
    var_threshold         = [.05 20]; % acceptable limits for variance in an epoch, relative to median of all epochs
    bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
    bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
    verbose               = true;
    
    [ts(:,:,data_channels), badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
        var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
    
    ts = ts(:,~badEpochs,:);
    ts(:,:, badChannels) = NaN;
    
    conditions = conditions(~badEpochs);
    
    
    %% Denoise data by regressing out nuissance channel time series
    
    % TODO: check whether this runs with NaNs in ts
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels
        if exist('./denoised_with_nuissance_data.mat', 'file')
            load(fullfile(data_pth{subject_num},'denoised_with_nuissance_data.mat'));
        else fprintf('Denoising data... \n');
            ts = meg_environmental_denoising(ts, environmental_channels,...
                data_channels, produce_figures);
        end
    end
    
    
    % --------------------------------------------------------------------
    % ------------------ ANALYZE THE PREPROCESSED DATA -------------------
    % --------------------------------------------------------------------
    %%
    % low pass filter the evoked signal
    %   Design a 70th order lowpass FIR filter with cutoff frequency of 75 Hz.
    df = designfilt('lowpassfir','FilterOrder',70,'CutoffFrequency',50, 'SampleRate', fs);
    D = mean(grpdelay(df)); % filter delay in samples
    
    ts_f = padarray(ts, [D 0 0], 0, 'post');
    ts_f = filter(df,ts_f);                   % Append D zeros to the input data
    ts_f = ts_f(D+1:end,:, :);                % Shift data to compensate for delay
    
    blank_trials     = conditions == find(blank_condition);
    evoked_grand     = squeeze(nanmean(ts_f(:, ~blank_trials, :) , 2));
    evoked_grand_sem = squeeze(nanstd(ts_f(:, ~blank_trials, :), [],  2)) / sqrt(sum(~blank_trials));
        
    evoked_grand_rms = sqrt(nansum(evoked_grand(:,data_channels).^2,2));
    [~, grand_idx] = max(evoked_grand_rms);
    peak_window = [-30:30] + grand_idx;
    figure; plot(evoked_grand_rms); hold on;
    % add plot lines for peak_window
    
    [~, which_time_point] = max(abs(evoked_grand(peak_window, data_channels)));
    which_time_point = which_time_point + min(peak_window) -1;
    
    linearInd = sub2ind([num_time_points size(evoked_grand,2)], which_time_point, data_channels);
    
    channel_peaks = abs(evoked_grand(linearInd)) ./ evoked_grand_sem(linearInd);
    
    figure;
    ft_plotOnMesh(channel_peaks); title('Evoked response (mean(peak)/sem(peak))');
    set(gca, 'CLim', [0 30]);
    
    
    figure;
    noise_pool = channel_peaks < 5;
    ft_plotOnMesh(double(noise_pool), [], [], [], 'interpolation', 'nearest'); title('noise pool');
                
end


