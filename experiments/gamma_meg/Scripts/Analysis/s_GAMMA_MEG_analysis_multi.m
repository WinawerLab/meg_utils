% s_GAMMA_MEG_analysis
%
% NOTE: A version of this script called HPC_Gamma runs on the HPC.
% If changes are made to this script, parallel changes should be made to
% HPC_Gamma
%
% After running this script or HPC_Gamma, results can be visualized with
% s_GAMMA_MEG_visualize.m
%
% Analyze data from MEG Gamma experiments. Subjects saw
% several kinds of stimuli, including gratings of various spatial
% frequencies, plaids, and noise patterns (phase-randomized with 1/f^n
% spectral power distrubutions)
%
% Stimuli were static, and on the screen either for 500 ms (subjects 1-3)
% or 1000 ms (subjects 4-6) with 500 ms ISI.
%
% Spectral data from each channel for each stimulus type are modeled as a
% mixture of a line and gaussian in log power / log frequency (meaning, a
% power law and a narrowband response)
%
% See also s_GAMMA_MEG_visualize, HPC_Gamma
%
% TODO:
%  MAYBE: source localize the signals

% Analysis options
%% Set analysis variables


data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

denoise_with_nonphys_channels = true;        % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

nboot                         = 1;           % number of bootstrap samples

produce_figures               = true;        % If you want figures in case of debugging, set to true

%denoise_via_pca               = false;       % Do you want to use meg denoise?
use_denoised                  = true;

fs                            = 1000;        % sample rate
epoch_start_end               = [0.050 1.049];% start and end of epoch, relative to trigger, in seconds

%intertrial_trigger_num        = 11;          % the MEG trigger value that corresponds to the intertrial interval
intertrial_trigger_num        = 14;          % the MEG trigger value that corresponds to the intertrial interval

save_images                   = false;
save_spectral_data            = false;

%which_sessions_to_analyze    = [16];    % subject 99 for synthetic data
which_sessions_to_analyze    = 18; % face / house / grating experiment

suffix                        = 'localregression_multi_10_boots_denoised';

%% Add paths
% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})


%% Loops over datasets
for session_num = which_sessions_to_analyze
    
    condition_names  = gamma_get_condition_names(session_num);
    
    blank_condition = strcmpi(condition_names, 'blank');
        
    if session_num == 99
        [ts, conditions] = gamma_synthetize_validation_data();
        
    else
        path_to_data = meg_gamma_get_path(session_num);

        save_pth = fullfile(path_to_data, 'processed');
        
        if ~exist(save_pth, 'dir'), mkdir(save_pth); end
        
        % --------------------------------------------------------------------
        % ------------------ PREPROCESS THE DATA -----------------------------
        % --------------------------------------------------------------------
        %% Load data (SLOW)
        if use_denoised
            denoised_data = load(fullfile(path_to_data, 'raw', sprintf('s0%d_denoisedData.mat',session_num)));
            tmp = denoised_data.denoisedts{1}; 
            raw_ts = permute(tmp,[2,3,1]);
            raw_ts = meg_load_sqd_data(fullfile(path_to_data, 'raw'), '*Gamma*');
        else
            raw_ts = meg_load_sqd_data(fullfile(path_to_data, 'raw'), '*Gamma*');
        end
        %% Extract triggers
        trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
        
        %% Make epochs
        [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
        
        clear raw_ts;
        % remove intertrial intervals
        iti               = conditions == intertrial_trigger_num;
        ts                = ts(:,~iti, :);
        conditions        = conditions(~iti);
    end
    % There are some weird unrelated triggers in the data, here we just
    % eliminate these.
    if sum(conditions == 15) > 0;
        idx           = find(conditions==15);
        ts(:,idx, :)  = [];
        conditions(idx) = [];
    end
    
    if sum(conditions == 12)  > 0 && length(condition_names) < 12;
        idx           = find(conditions==12);
        ts(:,idx, :)  = [];
        conditions(idx) = [];
    end
    
    
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
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels && session_num ~= 99
         fprintf('Denoising data... \n');
            ts = meg_environmental_denoising(ts, environmental_channels,...
                data_channels, produce_figures);        
    end
    
    
    % --------------------------------------------------------------------
    % ------------------ ANALYZE THE PREPROCESSED DATA -------------------
    % --------------------------------------------------------------------
    %% Spectral analysis
    
    % compute spectral data
    t = (1:size(ts,1))/fs;
    f = (0:length(t)-1)/max(t);
    
    spectral_data = abs(fft(ts))/length(t)*2;
    spectral_data_boots = zeros(size(ts,1), length(conditions_unique), length(data_channels), nboot);
    
    % compute the mean amplitude spectrum for each electrode in each condition
    fprintf('Computing bootstraps for each condition\n');
    for ii = 1:length(conditions_unique)
        fprintf('Condition %d of %d\n', ii, length(conditions_unique)); drawnow;
        
        % Binary vector to identify epochs with this condition
        these_epochs = conditions == conditions_unique(ii);
        
        % spectral data, time points x epochs x channel
        these_data = spectral_data(:,these_epochs,data_channels);
        
        if nboot > 1
            % reshape so that epochs are in rows (needed for bootstrp)
            these_data = permute(these_data, [2 1 3]);
            
            % log normalized mean of spectral power
            bootfun = @(x) squeeze(exp(nanmean(log(x),1)));
            
            % bootstat by definition is a matrix: nboot x (freq x channel)
            bootstat = bootstrp(nboot, bootfun, these_data);
            
            % reshape bootstat to 3D-array: nboot x freq x channel
            bootstat = reshape(bootstat, nboot, length(t), []);
            
            % spectral_data_boots is freq x condition x channel x boot
            spectral_data_boots(:,ii,:,:) = permute(bootstat,[2 3 1]);
            
        else
            spectral_data_boots(:,ii,:,:) = exp(nanmean(log(these_data),2));
        end
    end
    fprintf('Done!\n');
    
    % Summarize bootstrapped spectral by mean and std over bootstraps
    if save_spectral_data
        save(fullfile(save_pth,sprintf('spectral_data_%s.mat',suffix)),'spectral_data_boots')
    end
    
    spectral_data_mean = mean(spectral_data_boots, 4);
    
    %% Broadband and Gaussian Fit
    
    % Convert the amplitude spectrum in each channel and each epoch into 2
    % numbers, one for broadband and one for gamma
    
    %     f_use4fit = f((f>=35 & f < 40) |(f > 40 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
    
    % Note: smaller drop out of line noise frequency, in order to improve
    % modelfit
    f_use4fit = f((f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
    f_sel=ismember(f,f_use4fit);
    num_time_points = round((epoch_start_end(2)-epoch_start_end(1)+0.001)*fs);
    
    num_channels = length(data_channels);
%     fit_bl  = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
    w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
    w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
    gauss_f = NaN(num_channels, 1, nboot);                 % gaussian peak frequency
    fit_f2  = NaN(num_conditions,num_time_points,num_channels, nboot); % fitted spectrum
    
    warning off 'MATLAB:subsassigndimmismatch'
    
    % For each channel, fit each condition separatley
    fprintf('Fitting gamma and broadband values for each channel and each condition')
    % Fit each channel separately
    
    for chan = data_channels
        
        fprintf('Channel %d of %d\n', chan, length(data_channels)); drawnow;
        
        
        for bootnum = 1:nboot
            
            % the baseline is the same for all conditions, so just compute it once
            data_base = exp(mean(log(spectral_data_boots(:,:,chan, bootnum)),2));
            data_base = data_base';
            
            if all(isfinite(data_base))
                                
                data_fit = spectral_data_boots(:,:,chan, bootnum);
                
                [...
                    fit_bl(chan, :, bootnum), ...
                    w_pwr(chan, :, bootnum), ...
                    w_gauss(chan, :, bootnum),...
                    gauss_f(chan, 1, bootnum),...
                    fit_f2(:,:, chan, bootnum)] = ...
                    gamma_fit_data_localregression_multi(f,f_use4fit,data_base,data_fit);
                
%                                 fH = figure(1); clf;
%                                 plot(f, exp(fit_f2(:,:, chan, bootnum))','r', f, spectral_data_boots(:,:,chan, bootnum), 'k')
%                                 set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [10 200])
%                                 title(sprintf('Channel %d', chan))
%                                 pause(0.1)
            end
        end
    end
    fprintf('done!\n')
    
    warning on 'MATLAB:subsassigndimmismatch'
    
    gauss_f = repmat(gauss_f, [1 num_conditions 1]);
    
    
    % summarize bootstrapped fits
    fit_bl_mn  = nanmean(fit_bl,3);
    w_pwr_mn   = nanmean(w_pwr,3);
    w_gauss_mn = nanmean(w_gauss,3);
    gauss_f_mn = nanmean(gauss_f,2);
    fit_f2_mn  = nanmean(fit_f2,4);
    
    
    %% Save Processed Data
    filename = fullfile(save_pth, sprintf('s0%d_%s.mat',session_num,suffix));
    save (filename, 'num_conditions', 'f_sel', 'data_channels', 'nboot', 'f_use4fit', ...
        'fit_bl', 'w_pwr', 'w_gauss', 'gauss_f', 'fit_f2', 'w_gauss_mn', 'w_pwr_mn', 'fit_bl_mn');
    
    
    
end


