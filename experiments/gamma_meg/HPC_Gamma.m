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
for subject_num = which_data_sets_to_analyze
    
    if subject_num == 5 || subject_num == 6 
        intertrial_trigger_num = 11; % the MEG trigger value that corresponds to the intertrial interval
        epoch_start_end        = [0.55 1.049];% start and end of epoch, relative to trigger, in seconds

    else
        intertrial_trigger_num = 10;
        epoch_start_end        = [0.05 .549];% start and end of epoch, relative to trigger, in seconds

    end
    
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
    if intertrial_trigger_num == 11 % in this case, the blanks in between images have trigger number 11
        iti               = conditions == intertrial_trigger_num;
        ts                = ts(:,~iti, :);
        conditions        = conditions(~iti);
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
        verbose               = true;
        
        [ts(:,:,data_channels), badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
            var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
        
        ts = ts(:,~badEpochs,:);
        ts(:,:, badChannels) = NaN;
        
        conditions = conditions(~badEpochs);
        
    end
    
    
    %% Denoise data by regressing out nuissance channel time series
    
    % TODO: check whether this runs with NaNs in ts
    
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
    spectral_data_mean = mean(spectral_data_boots, 4);
    spectral_data_std  =  std(spectral_data_boots, [], 4);
    spectral_data_snr   = spectral_data_mean./spectral_data_std;
    
    %% Broadband and Gaussian Fit
    
    % Convert the amplitude spectrum in each channel and each epoch into 2
    % numbers, one for broadband and one for gamma
    
    f_use4fit = f((f>=35 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
    f_sel=ismember(f,f_use4fit);
%     num_time_points = double(epoch_start_end(2)-epoch_start_end(1))*1000;

    
    num_channels = length(data_channels);
    out_exp = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
    w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
    w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
    gauss_f = NaN(num_channels,num_conditions, nboot);     % gaussian peak frequency
    fit_f2  = NaN(num_conditions,500,num_channels, nboot); % fitted spectrum
    
    warning off 'MATLAB:subsassigndimmismatch'
    
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
    warning on 'MATLAB:subsassigndimmismatch'
    
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
    parsave([fname '.mat'], 'out_exp', out_exp, 'w_pwr', w_pwr, ...
        'w_gauss', w_gauss, 'gauss_f', gauss_f,...
        'fit_f2', fit_f2, 'nboot', nboot);

    
    
end

function ts = meg_remove_bad_epochs(outliers, ts)
% epochs x channel
num_time_points = size(ts,1);
num_epochs      = size(ts,2);
num_channels    = size(ts,3);

ts = reshape(ts, [num_time_points, num_epochs*num_channels]);

ts(:, logical(outliers(:))) = NaN;

ts = reshape(ts, [num_time_points, num_epochs, num_channels]);

return


