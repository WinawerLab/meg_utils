function HPC_Gamma(which_data_sets_to_analyze, nboot)

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
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

produce_figures               = true;        % If you want figures in case of debugging, set to true

denoise_via_pca               = false;       % Do you want to use megdenoise?

fs                            = 1000;        % sample rate
epoch_start_end               = [0.550 1.05];% start and end of epoch, relative to trigger, in seconds

intertrial_trigger_num        = 11;          % the MEG trigger value that corresponds to the intertrial interval

save_images                   = false;
verbose                       = false;

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
parfor subject_num = which_data_sets_to_analyze
    
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
    
    %% Find bad epochs
    if remove_bad_epochs
        
        % This identifies any epochs whos variance is outside some multiple of the
        % grand variance
        bad_epochs = meg_find_bad_epochs(ts(:,:,data_channels), [.05 20]);
        
        % any epoch in which more than 10% of channels were bad should be removed
        % entirely
        epochs_to_remove = mean(bad_epochs,2)>.1;
        
        % once we remove 'epochs_to_remove', check whether any channels have more
        % than 10% bad epochs, and we will remove these
        channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.1;
        
        bad_epochs(epochs_to_remove,:) = 1;
        bad_epochs(:,channels_to_remove) = 1;
        
        if verbose
        figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')
        end
        
        ts = meg_remove_bad_epochs(bad_epochs, ts);
    end
    
    
    %% Denoise data by regressing out nuissance channel time series
    
    % TODO: check whether this runs with NaNs in ts
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels
        if exist('./denoised_with_nuissance_data.mat', 'file')
            load(fullfile(data_pth{subject_num},'denoised_with_nuissance_data.mat'));
        else
            if verbose; fprintf('Environmentally denoise data.. This may take a couple of seconds\n'); end
            ts = meg_environmental_denoising(ts, environmental_channels,...
                data_channels, false);
        end
    end
    
    
    % --------------------------------------------------------------------
    % ------------------ ANALYZE THE PREPROCESSED DATA -------------------
    % --------------------------------------------------------------------
    %% Spectral analysis
    
    % compute spectral data
    t = (1:size(ts,1))/fs;
    f = (0:length(t)-1)/max(t);
    nboot = 50; % number of bootstrap samples
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
    spectral_data_mean = mean(spectral_data_boots, 4);
    spectral_data_std  =  std(spectral_data_boots, [], 4);
    spectral_data_snr   = spectral_data_mean./spectral_data_std;
    
    %% Broadband and Gaussian Fit
    
    % Convert the amplitude spectrum in each channel and each epoch into 2
    % numbers, one for broadband and one for gamma
    
    f_use4fit = f((f>=35 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
    f_sel=ismember(f,f_use4fit);
    
    num_channels = length(data_channels);
    out_exp = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
    w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
    w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
    gauss_f = NaN(num_channels,num_conditions, nboot);     % gaussian peak frequency
    fit_f2  = NaN(num_conditions,500,num_channels, nboot); % fitted spectrum
    
%     warning off 'MATLAB:subsassigndimmismatch'
    warning off
    
    % For each channel, fit each condition separatley
    if verbose; fprintf('Fitting gamma and broadband values for each channel and each condition'); end;
    for cond = 1:num_conditions
        if verbose; fprintf('.'); drawnow; end;
        % Fit each channel separately
        for chan = data_channels
            
            
            for bootnum = 1:nboot
                
                data_fit  = spectral_data_boots(:,cond,chan, bootnum);
                data_base = spectral_data_boots(:,blank_condition,chan, bootnum);
                
                % try/catch because bad channels / bad epochs were replaced by
                % NaNs, and NaNs will cause an error
                try
                    [...
                        out_exp(chan, cond, bootnum), ...
                        w_pwr(chan, cond, bootnum), ...
                        w_gauss(chan, cond, bootnum),...
                        gauss_f(chan, cond, bootnum),...
                        fit_f2(cond,:, chan, bootnum)] = ...
                        gamma_fit_data(f,f_use4fit,data_base,data_fit);
                catch ME
                    warning(ME.identifier, ME.message)
                end
            end
        end
    end
    if verbose; fprintf('done!\n'); end;
    
%     warning on 'MATLAB:subsassigndimmismatch'
    warning on
    
    % summarize bootstrapped fits
    out_exp_mn = nanmean(out_exp,3);
    w_pwr_mn   = nanmean(w_pwr,3);
    w_gauss_mn = nanmean(w_gauss,3);
    gauss_f_mn = nanmean(gauss_f,3);
    fit_f2_mn  = nanmean(fit_f2,4);
    
    out_exp_sd = nanstd(out_exp,[],3);
    w_pwr_sd   = nanstd(w_pwr,[],3);
    w_gauss_sd = nanstd(w_gauss,[],3);
    gauss_f_sd = nanstd(gauss_f,[],3);
    fit_f2_sd  = nanstd(fit_f2,[],4);
    
    % Save data
    fname = sprintf(fullfile(rootPath,'HPC','Data','s0%d_bootstrappedData'),which_data_sets_to_analyze);
    parsave([fname '.mat'], 'out_exp_mn', out_exp_mn, 'w_pwr_mn', w_pwr_mn, ...
        'w_gauss_mn', w_gauss_mn, 'gauss_f_mn', gauss_f_mn,...
        'fit_f2_mn', fit_f2_mn, ...
        'out_exp_sd', out_exp_sd, 'w_pwr_sd', w_pwr_sd, ...
        'w_gauss_sd', w_gauss_sd, 'gauss_f_sd', gauss_f_sd,...
        'fit_f2_sd', fit_f2_sd, 'nboot', nboot);

    
    
end

return


