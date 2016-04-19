function results = gamma_spectral_analysis(ts, conditions, nboot, sessionNum)
%% [spectralData, fitData, summaryStats] = gamma_spectral_analysis(ts, conditions)
%  performs the spectral analysis, bootstrapping, and curve fitting for
%  s_GAMMA_pipeline
% input: ts - epoched timeseries, can be denoised or not yet denoised
%      : conditions - a vector of integers representing the condition of
%      each epoch
% output: spectralResults - a structure containing:
%            1) 'f_sel'
%            2) 'f_use4fit'
%            3) 'fit_bl'
%            4) 'w_pwr'
%            5) 'w_gauss'
%            6) 'gauss_f'
%            7) 'fit_f2'
%            8) 'spectral_data'

%% Parameters
data_channels = 1:157;
fs = 1000;
epoch_start_end = [0.050 1.049];

save_spectral_data = false; % all the bootstraps 
save_results       = true; % spectral mean and fits
suffix = datestr(now, 'mm.dd.yy');

%% Calculate Spectral Data



% remove ITI epochs from data
ITInum = max(conditions);
ts = ts(:,conditions ~= ITInum, :);
ITI = conditions == ITInum;
conditions = conditions(~ITI);

conditions_unique = unique(conditions);
num_conditions = length(conditions_unique);

% compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);

% use these frequencies
fitFreq = f((f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));

spectral_data = abs(fft(ts))/length(t)*2;
% freq x conditions x channels x bootstraps
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
spectral_data_mean = nanmean(spectral_data_boots, 4);


%% Broadband and Gaussian Fit

% binary array conrresponding to freq
f_sel = ismember(f,fitFreq);
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
                gamma_fit_data_localregression_multi(f,fitFreq,data_base,data_fit);
            
            %fH = figure(1); clf;
            %plot(f, exp(fit_f2(:,:, chan, bootnum))','r', f, spectral_data_boots(:,:,chan, bootnum), 'k')
            %set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [10 200])
            %title(sprintf('Channel %d', chan))
            %pause(0.1)
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

results = struct('spectral_data_mean', spectral_data_mean, 'fit_bl_mn', fit_bl_mn, 'w_pwr_mn', w_pwr_mn, 'w_gauss_mn', w_gauss_mn, ...
    'gauss_f_mn', gauss_f_mn, 'fit_f2_mn', fit_f2_mn);

save_pth = fullfile(meg_gamma_get_path(sessionNum), 'processed');

if save_spectral_data
    save(fullfile(save_pth,sprintf('spectral_data_%s.mat',suffix)),'spectral_data_boots', '-v7.3')
end

if save_results
    save(fullfile(save_pth,sprintf('s%03d_summary&fits_%s.mat', sessionNum, suffix)),'results', '-v7.3')
end


end