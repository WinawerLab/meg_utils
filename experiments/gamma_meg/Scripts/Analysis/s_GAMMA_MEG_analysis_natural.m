%% OBSOLETE

error('This script is obsolete')

return

%% s_GAMMA_MEG_analysis_natural.m
% Analyse data from MEG gamma experiments
% (Written with the natural stimuli data in mind, but should be
% multipurpose)
% Obtains spectral data: 
%
% TODO:
%   - create 'condition' file specifying each experiment's conditions
%   - analysis
%   - denoise
%   - plot

%% Set analysis variables
% paths
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                   = '*_Gamma_*subj*';
meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

% channel assignment
data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

% experimental parameters
fs                            = 1000; %1000 Hz sampling rate
epoch_start_end               = [0.050 1.049]; %start/end of epoch relative to trigger in seconds
intertrial_trigger_num        = 14; %trigger number of intertrial interval

% processing options
denoise_with_nonphys_channels = true;        % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range
denoise_via_pca               = false;

% analysis options
nboot                         = 1;
produce_figures               = true;
save_images                   = false;
save_spectral_data            = true;

which_data                    = 19;

suffix                        = 'preliminary';

%% Paths
d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

%% loop over datasets
for session_num = which_data
    %% get conditions
    condition_names  = gamma_get_condition_names(session_num-1);
    blank_condition = strcmpi(condition_names, 'blank');
    
    % define and make directory for saved data/images
    save_pth = fullfile(project_pth, 'Images', subj_pths{session_num-1});
    if ~exist(save_pth, 'dir'), mkdir(save_pth); end
    
    %% Load (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{session_num-1}, 'raw'), '*Gamma*');
    
    %% Extract Trigger Sequence
    trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    
    %% Make Epochs
    [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
    % remove intertrial intervals
    iti               = conditions == intertrial_trigger_num;
    ts                = ts(:,~iti, :);
    conditions        = conditions(~iti);
    
    conditions_unique = unique(conditions);
    num_conditions    = length(condition_names);
    
    %% Remove Bad Epochs
    var_threshold    = [0.05 20];
    bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
    bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
    verbose               = true;
    
    [ts(:,:,data_channels), badChannels, badEpochs] = meg_preprocess_data(ts(:,:,data_channels), ...
        var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
    ts = ts(:,~badEpochs,:);
    ts(:,:, badChannels) = NaN; % turn data in bad epoch/channel = NaNs
    
    conditions = conditions(~badEpochs);
    
     %% Spectral analysis
    
    % compute spectral data
    t = (1:size(ts,1))/fs;
    f = (0:length(t)-1)/max(t);
    
    spectral_data = abs(fft(ts))/length(t)*2; % f x epochs x channel
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
        save(fullfile(project_pth, subj_pths{session_num-1},'processed',sprintf('spectral_data_%s.mat',suffix)),'spectral_data_boots')
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
    out_exp = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
    w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
    w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
    gauss_f = NaN(num_channels,num_conditions, nboot);     % gaussian peak frequency
    fit_f2  = NaN(num_conditions,num_time_points,num_channels, nboot); % fitted spectrum
    
    warning off 'MATLAB:subsassigndimmismatch'
    
    % For each channel, fit each condition separatley
    fprintf('Fitting gamma and broadband values for each channel and each condition')
    % Fit each channel separately
    
    fH = figure(1); clf; 
    for chan = data_channels
        
        fprintf('Channel %d of %d\n', chan, length(data_channels)); drawnow;
        
        
        for bootnum = 1:nboot
            
            % the baseline is the same for all conditions, so just compute it once
            data_base = exp(mean(log(spectral_data_boots(:,:,chan, bootnum)),2));
            data_base = data_base';
            
            if all(isfinite(data_base))
                
                for cond = 1:num_conditions
                    
                    data_fit = spectral_data_boots(:,cond,chan, bootnum);
                    data_fit = data_fit';
                    
                    % try/catch because bad channels / bad epochs were replaced by
                    % NaNs, and NaNs will cause an error
                    
                    [...
                        ~, ...
                        w_pwr(chan, cond, bootnum), ...
                        w_gauss(chan, cond, bootnum),...
                        gauss_f(chan, cond, bootnum),...
                        fit_f2(cond,:, chan, bootnum)] = ...
                        ...      gamma_fit_data(f,f_use4fit,data_base,data_fit);
                        gamma_fit_data_localregression(f,f_use4fit,data_base,data_fit);
                end
                
                plot(f, exp(fit_f2(:,:, chan, bootnum))','r', f, spectral_data_boots(:,:,chan, bootnum), 'k')
                set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [10 200])
                title(sprintf('Channel %d', chan))
                drawnow
            end
        end
    end
    fprintf('done!\n')
    
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
    
    out_exp_md = nanmedian(out_exp,3);
    w_pwr_md   = nanmedian(w_pwr,3);
    w_gauss_md = nanmedian(w_gauss,3);
    gauss_f_md = nanmedian(gauss_f,3);
    fit_f2_md  = nanmedian(fit_f2,4);
    
    %% Save Processed Data
    filename = fullfile(project_pth, subj_pths{session_num-1}, 'processed', sprintf('s0%d_%s.mat',session_num,suffix));
    save (filename, 'project_pth', 'num_conditions', 'f_sel', 'data_channels', 'nboot', 'f_use4fit', ...
        'out_exp', 'w_pwr', 'w_gauss', 'gauss_f', 'fit_f2', 'w_gauss_mn', 'w_pwr_mn');
end

