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
colormap(params.baselineCondition,:) = [0 0 0];

% Define the prefix and suffix of files to be saved
if SAVE_FIGS
    isDenoised = '';
    if param.pcaDenoise, isDenoised = '_denoised'; end
    prefix = sprintf('s_%3d', sessionNum);
    suffix = sprintf('_%dbootstraps%s_%s', nBoot, isDenoised,...
        datestr(now, 'mm_dd_yy'));
end
%% 1. All conditions plotted on the same spectrogram for each channel
f_plot = f;
f_plot(~f_sel) = NaN;

% fit_f2_mn is cond x freq x chan
fit_f2_mn = exp(permute(results.fit_f2_mn, [2 1 3])); % NOTE: spectral data is exponentiated in gamma_spectral_analysis

fig3 = figure; set(gcf, 'color', 'w');

plotThis = true;
if plotThis
    for chan = 1:157
        subplot(1,2,1); cla
        set(gca, 'Colororder', colormap, 'XScale', 'log', 'XLim', [35 200]); hold on
        plot(f_plot, specData(:, :, chan),...
            'LineWidth', 2); %'color', colormap(ii,:,:)
        legend(conditionNames)
        yl = get(gca, 'YLim');
        
        subplot(1,2,2); cla
        set(gca, 'Colororder', colormap, 'XScale', 'log', 'XLim', [35 200]); hold on
        plot(f_plot, fit_f2_mn(:, :, chan),...
            'LineWidth', 2); % 'color', colormap(ii,:,:),
        set(gca, 'YLim', yl);
        waitforbuttonpress;
    end
end
%% 2. All conditions' line fit plotted on the same spectrogram for each channel

% fit_f2_mn is chan x freq x chan
% fit_f2_mn = permute(results.fit_f2_mn, [2 1 3]);
% 
% 
% for chan = 1:157
%     figure(1); clf;
%     set(gcf, 'color', 'w');
%     hold all;
%     
%     for ii = 1:length(conditionNames)
%         plot(f_plot, smooth(fit_f2_mn(:, ii, chan), 2)',...
%             'color', colormap(ii,:,:),'LineWidth', 2);
%     end
%     legend(conditionNames)
%     title(chan)
%     waitforbuttonpress;
% end

%% All channels, fitted spectrogram for each condition

clf;
figure(3);
for chan = 1:10
    set(gcf, 'Name', sprintf('Channel %d', chan));
   for cond = 1:length(conditionNames)-1
       subplot(4,3,cond); cla;
       % plot given stimuli condition
       plot(f_plot, fit_f2_mn(:,cond,chan), 'Color', colormap(cond,:)); hold on;
       % plot baseline fit
       plot(f_plot, fit_f2_mn(:,length(conditionNames), chan), 'Color', 'k');
       hold off;
       %set(gca, 'Color', colormap(cond, :));
       title(cell2mat(conditionNames(cond)));
       %title(sprintf('%s BB:%f.3;  G:%f.3', cell2mat(conditionNames(cond)), results.w_pwr_mn(chan, cond), results.w_gauss_mn(chan, cond)));
       set(gca, 'XScale', 'log', 'XLim', [35 200]);
   end
   waitforbuttonpress;


end

