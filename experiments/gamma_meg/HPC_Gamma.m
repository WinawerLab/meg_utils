function HPC_Gamma(which_sessions_to_analyze, nboot)
%
% NOTE:  A version of this function called s_GAMMA_MEG_analysis runs
% without the HPC. If changes are made to this script, parallel changes
% should be made to s_GAMMA_MEG_analysis

% Function to analyze Gamma analysis on the HPC

% Analysis options
%% Set analysis variables
rootPath                      = which('HPC_Gamma');
rootPath                      = fileparts(rootPath);

project_pth                   = fullfile(rootPath,'HPC','Data');

% data to be analysed
data_pth                      = '*_Gamma_*subj*';

data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

denoise_with_nonphys_channels = true;        % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold

fs                            = 1000;        % sample rate

save_spectral_data            = true;

suffix                        = 'localregression_multi';

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

blank_condition = strcmpi(condition_names, 'blank');
%% Add paths

%change server-1 back to server
% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip', 'yokogawa_defaults')

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);
%% Loops over datasets
for session_num = which_sessions_to_analyze
    verbose                       = true;
%     if session_num > 4
        intertrial_trigger_num = 11; % the MEG trigger value that corresponds to the intertrial interval
        epoch_start_end        = [0.05 1.049]; %[0.55 1.049];% start and end of epoch, relative to trigger, in seconds
        
%     else
%         intertrial_trigger_num = 10;
%         epoch_start_end        = [0.05 .549];% start and end of epoch, relative to trigger, in seconds
%         
%     end
   
    
    % --------------------------------------------------------------------
    % ------------------ PREPROCESS THE DATA -----------------------------
    % --------------------------------------------------------------------
    %% Load data (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{session_num}, 'raw'), '*Gamma*');
    
    %% Extract triggers
    trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    
    %% Make epochs
    [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
    % remove intertrial intervals
    if intertrial_trigger_num == 11 % in this case, the blanks in between images have trigger number 11
        iti               = conditions == intertrial_trigger_num;
        ts                = ts(:,~iti, :);
        conditions        = conditions(~iti);
        
        % There are some weird unrelated triggers in the data, here we just
        % eliminate these.
        if sum(conditions == 15) > 0;
            idx           = find(conditions==15);
            ts(:,idx, :)  = [];
            conditions(idx) = [];
        end
        
        if sum(conditions == 12) > 0;
            idx           = find(conditions==12);
            ts(:,idx, :)  = [];
            conditions(idx) = [];
        end
        
    else % in this case, the blanks in between images have trigger number 10
        ts                = ts(:,1:2:end, :);
        conditions        = conditions(1:2:end);
    end
    
    
    conditions_unique = unique(conditions);
    num_conditions    = length(condition_names);
    

    
    
    %% Find bad epochs
    if remove_bad_epochs
        
        var_threshold         = [.05 20]; % acceptable limits for variance in an epoch, relative to median of all epochs
        bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
        bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
        verbose               = false;
        
        [ts(:,:,data_channels), badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
            var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
        
        ts = ts(:,~badEpochs,:);
        ts(:,:, badChannels) = NaN;
        
        conditions = conditions(~badEpochs);
        
    end
    
        %% Denoise data by regressing out nuissance channel time series
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels
        if verbose; fprintf('Environmentally denoise data.. This may take a couple of seconds\n'); end
        ts = meg_environmental_denoising(ts, environmental_channels,...
            data_channels, false);
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
    if verbose; fprintf('Computing bootstraps for each condition'); end;
    for ii = 1:length(conditions_unique)
        if verbose; fprintf('.'); drawnow; end;
        
        % Binary vector to identify epochs with this condition
        these_epochs = conditions == conditions_unique(ii);
        
        % spectral data, time points x epochs x channel
        these_data = spectral_data(:,these_epochs,data_channels);
        
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
        
    end
    if verbose; fprintf('Done!\n'); end;
    
    
    % Summarize bootstrapped spectral by mean and std over bootstraps
    if save_spectral_data
%         if ~exist(fullfile(project_pth, subj_pths{session_num},'processed'),'dir'); mdkir(fullfile(project_pth, subj_pths{session_num},'processed'));end;
        parsave(fullfile(project_pth, subj_pths{session_num},sprintf('spectral_data_%s.mat',suffix)),'spectral_data_boots', spectral_data_boots)
    end
    
    spectral_data_mean = mean(spectral_data_boots, 4);
    %% Broadband and Gaussian Fit
    
    % Convert the amplitude spectrum in each channel and each epoch into 2
    % numbers, one for broadband and one for gamma
    
    %     f_use4fit = f((f>=35 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
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
                    gauss_f(chan, :, bootnum),...
                    fit_f2(:,:, chan, bootnum)] = ...
                    gamma_fit_data_localregression_multi(f,f_use4fit,data_base,data_fit);
                
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
    
    % Save data
    fname = fullfile(rootPath,'HPC','Data',sprintf('s0%d_%s',session_num,suffix));
    parsave([fname '.mat'], 'fit_bl', fit_bl, 'w_pwr', w_pwr, ...
        'w_gauss', w_gauss, 'gauss_f', gauss_f,...
        'w_gauss_mn', w_gauss_mn, 'w_pwr_mn', w_pwr_mn,  'fit_bl_mn', fit_bl_mn, ...
        'fit_f2', fit_f2, 'nboot', nboot, 'f_use4fit', f_use4fit, 'f_sel',f_sel);
   
    
    
end


