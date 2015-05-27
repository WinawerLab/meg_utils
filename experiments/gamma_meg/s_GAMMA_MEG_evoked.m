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
    'Plaid'...
    'Blank'};

which_data_sets_to_analyze = 5;
blank_condition = strcmpi(condition_names, 'blank');

%% Add paths

%change server-1 back to server
meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

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
    
    evoked = zeros(size(ts,1), num_conditions, size(ts,3));
    evoked_sem = zeros(size(ts,1), num_conditions, size(ts,3));

    for c = 1:length(conditions_unique)
       evoked(:,c, :) = nanmean(ts_f(:, conditions == c, :) , 2); 
       evoked_sem(:,c, :) = nanstd(ts_f(:, conditions == c, :), [],  2) / sqrt(sum(conditions==c));
    end
    
    channel_to_plot = 47;
    %channel_to_plot = 133;
    figure(100), clf, hold all
    plot(bsxfun(@plus, evoked(:,:, channel_to_plot), 150*(1:num_conditions)), 'r')    
    plot(bsxfun(@plus, evoked(:,:, channel_to_plot)+evoked_sem(:,:, channel_to_plot), 150*(1:num_conditions)), 'k--')
    plot(bsxfun(@plus, evoked(:,:, channel_to_plot)-evoked_sem(:,:, channel_to_plot), 150*(1:num_conditions)), 'k--')
    set(gca, 'YTick', 150*(1:num_conditions), 'YGrid', 'on');
    
    
    
    %%
end


