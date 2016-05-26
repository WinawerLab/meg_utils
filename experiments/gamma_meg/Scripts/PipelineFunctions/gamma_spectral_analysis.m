function results = gamma_spectral_analysis(ts, params)
%% [spectralData, fitData, summaryStats] = gamma_spectral_analysis(ts, conditions)
%  performs the spectral analysis, bootstrapping, and curve fitting for
%  s_GAMMA_pipeline
% input: ts - epoched timeseries, can be denoised or not yet denoised
%      : params - structure containing experimental parameters obtained
%        from gamma_get_parameters.m
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
data_channels     = 1:157;
fs                = 1000;
epoch_start_end   = params.epochRange;
nBoot             = params.nBoot;
conditions        = params.conditionVector;
sessionNum        = params.sessionNumber;

% save options
SAVE_SPECTRAL_DATA = false; % all the bootstraps 
SAVE_RESULTS       = true; % spectral mean and fits

%% Get suffix for this spectral analysis based on parameters

thisDate = datestr(now, 'mm.dd.yy');

if params.pcaDenoise
    denoiseStr = '_denoised';
else
    denoiseStr = '';
end

suffix   = sprintf('%dboots%s_%s', nBoot, denoiseStr, thisDate);

%% Calculate Spectral Data
% remove ITI epochs from data

% ITI epochs already removed somewhere??
%  ITInum = params.ITI;
%  ts = ts(:,conditions ~= ITInum, :);
%  ITI = conditions == ITInum;
%  conditions = conditions(conditions~=ITI);


conditions_unique = unique(conditions);
num_conditions = length(conditions_unique);

% compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);

% use these frequencies
fitFreq = f((f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));

spectral_data = abs(fft(ts))/length(t)*2;
% freq x conditions x channels x bootstraps
spectral_data_boots = zeros(size(ts,1), length(conditions_unique), length(data_channels), nBoot);

% compute the mean amplitude spectrum for each electrode in each condition
fprintf('Computing bootstraps for each condition\n');
for ii = 1:length(conditions_unique)
    fprintf('Condition %d of %d\n', ii, length(conditions_unique)); drawnow;
    
    % Binary vector to identify epochs with this condition
    these_epochs = conditions == conditions_unique(ii);
    
    % spectral data, time points x epochs x channel
    these_data = spectral_data(:,these_epochs,data_channels);
    
    if nBoot > 1
        % reshape so that epochs are in rows (needed for bootstrp)
        these_data = permute(these_data, [2 1 3]);
        
        % log normalized mean of spectral power
        bootfun = @(x) squeeze(exp(nanmean(log(x),1)));
        
        % bootstat by definition is a matrix: nboot x (freq x channel)
        bootstat = bootstrp(nBoot, bootfun, these_data);
        
        % reshape bootstat to 3D-array: nboot x freq x channel
        bootstat = reshape(bootstat, nBoot, length(t), []);
        
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

% binary array describing the frequiencies not ommitted from the fit
params.f_sel = ismember(f,fitFreq);
num_time_points = round((epoch_start_end(2)-epoch_start_end(1)+0.001)*fs);

% per loop sizes
% fit_bl: chan x 112 x boot
% w_pwr: chan x cond x boot
% gauss_f : chan x 1 x boot
% fit_f2: freq x cond x chan x boot

num_channels = length(data_channels);
% fit_bl  = NaN(num_channels,length(fitFreq), nBoot);     % slope of spectrum in log/log space
w_pwr   = NaN(num_channels,num_conditions, nBoot);     % broadband power
w_gauss = NaN(num_channels,num_conditions, nBoot);     % gaussian height
gauss_f = NaN(num_channels, 1, nBoot);                 % gaussian peak frequency
fit_f2  = NaN(num_conditions,800,num_channels, nBoot); % fitted spectrum

warning off 'MATLAB:subsassigndimmismatch'

% For each channel, fit each condition separatley
fprintf('Fitting gamma and broadband values for each channel and each condition')
% Fit each channel separately

for chan = data_channels
    
    fprintf('Channel %d of %d\n', chan, length(data_channels)); drawnow;
        
    for bootnum = 1:nBoot
        
        % the baseline is the same for all conditions, so just compute it once
        data_base = exp(mean(log(spectral_data_boots(:,:,chan, bootnum)),2));
        data_base = data_base';
        
        if all(isfinite(data_base)) && mean(data_base) > 0
            
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
    'gauss_f_mn', gauss_f_mn, 'fit_f2_mn', fit_f2_mn, 'fitFreq', fitFreq, 'params', params);

save_pth = fullfile(meg_gamma_get_path(sessionNum), 'processed');

%% Save spectral data and results
% this is a really big file
if SAVE_SPECTRAL_DATA
    save(fullfile(save_pth,sprintf('s%03d_spectral_data_%s.mat', sessionNum, ...
        suffix)),'spectral_data_boots', '-v7.3')
end
% this is not
if SAVE_RESULTS
    save(fullfile(save_pth,sprintf('s%03d_summary&fits_%s.mat', sessionNum, ...
        suffix)),'results', '-v7.3')
end


end