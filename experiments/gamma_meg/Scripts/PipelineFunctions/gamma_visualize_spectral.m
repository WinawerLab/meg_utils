function [] = gamma_visualize_spectral(params, results)
% creates figures out of the spectral data obtained from
% gamma_spectral_analysis.m
%
% Input: params, a structure containing experimental and analysis parameters
%                such as condition names, epoch length, bootstraps, and how
%                the data is denoised. (See: gamma_get_parameters)
%        results, a structure containing: - spectral_data_mean (of each
%                                           condition and channel)
%                                         - fit_bl_mn (median baseline fit)
%                                         - w_pwr_mn (mean broadband
%                                           height of each cond x chan)
%                                         - w_gauss_mn (mean gamma height)
%                                         - gauss_f_mn (gamma peak freq)
%                                         - fit_f2_mn (fitted spectrum)
%
% Nicholas Chua 2016

%% Get parameters and spectral data

%options
SAVE_FIGS = false;

% get experiment params
sessionNum     = params.sessionNumber;
nBoot          = params.nBoot;
conditionNames = gamma_get_condition_names(sessionNum);

% get analysis results
specData       = results.spectral_data_mean;
numFreq        = size(specData, 1); % the number of frequency bins in data
%fitFreq = f((f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
fitFreq        = results.fitFreq; % the frequencies used to fit data

fs             = 1000; % resolution
t              = (1:numFreq)/fs;
f              = (0:length(t)-1)/max(t);
f_sel          = ismember(f, fitFreq); % boolean array of frequencies used


% get figure params
colormap = parula(length(conditionNames));

% Define the prefix and suffix of files to be saved
if SAVE_FIGS
    isDenoised = '';
    if param.pcaDenoise, isDenoised = '_denoised'; end
    prefix = sprintf('s_%3d', sessionNum);
    suffix = sprintf('_%dbootstraps%s_%s', nBoot, isDenoised,...
        datestr(now, 'mm_dd_yy'));
end
%% 1. All conditions plotted on the same spectrogram for each channel
plotThis = false;
if plotThis
    for chan = 1:157
        figure(1); clf;
        set(gcf, 'color', 'w');
        hold all;
        
        for ii = 1:length(conditionNames)
            plot(f(f_sel), smooth(specData(f_sel, ii, chan), 2)',...
                'color', colormap(ii,:,:),'LineWidth', 2);
        end
        legend(conditionNames)
        waitforbuttonpress;
    end
end
%% 2. All conditions' line fit plotted on the same spectrogram for each channel

% fit_f2_mn is chan x freq x chan
fit_f2_mn = permute(results.fit_f2_mn, [2 1 3]);


for chan = 1:157
    figure(1); clf;
    set(gcf, 'color', 'w');
    hold all;
    
    for ii = 1:length(conditionNames)
        plot(f(f_sel), smooth(fit_f2_mn(f_sel, ii, chan), 2)',...
            'color', colormap(ii,:,:),'LineWidth', 2);
    end
    legend(conditionNames)
    waitforbuttonpress;
end


end

