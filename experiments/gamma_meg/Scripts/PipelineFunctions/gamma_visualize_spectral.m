function [] = gamma_visualize_spectral(results, opt)
% creates figures out of the spectral data obtained from
% gamma_spectral_analysis.m
%
% Input: 
%   results (a struct which includes):
%       spectralDataMean    :    3D array containing spectral data averaged across bootstraps (time x conditions x channels) 
%       fitBaselineMean     :    2D array containing baseline fit averaged across bootstraps (channels x selected frequencies) 
%       broadbandPowerMean  :    2D array containing broadband power elevation fit averaged across bootstraps (channels x conditions) 
%       gammaPowerMean      :    2D array containing narrowband gamma power gaussian fit averaged across bootstraps (channels x conditions) 
%       gammaPeakFreqMean   :    1D array containing peak frequency of gaussian bump (channels x 1) 
%       modelFitAllFreqMean :    3D array containing modelfit for each frequency, condition, channel (frequencies, conditions, channels)
%       fitFreq             :    1D array containing all frequencies used for model fit
%       opt                 :    struct with options used when analyzing this dataset
%
% First version: Nicholas Chua April 2016
%       7.7.2016: Clean up (EK)

%% Prepare data
conditionNames = gamma_get_condition_names(opt.params.sessionNumber);
specData       = results.spectralDataMean;
numFreq        = size(specData, 1); % the number of frequency bins in data
t              = (1:numFreq)/opt.fs;
f              = (0:length(t)-1)/max(t);
f_sel          = ismember(f, opt.fitFreq); % boolean array of frequencies used

% Define colors
colormap = parula(length(conditionNames));
colormap(opt.params.baselineCondition,:) = [0 0 0];

% Define the file name of saved figures
if opt.saveFigures
    if opt.MEGDenoise; postFix = '_denoised'; else postFix = ''; end
    preFix = sprintf('s_%3d', opt.params.sessionNumber);
    saveName = sprintf('%s_boots%d%s_%s', preFix, opt.nBoot, postFix,...
        datestr(now, 'mm_dd_yy'));
end

%% 1. All conditions plotted on the same spectrogram for each channel

% Define frequencies to plot
f_plot = f;
f_plot(~f_sel) = NaN;

% Exponentiate modelfit modelFitAllFreqMean (conditions x frequencies x channels)
fit = exp(permute(results.modelFitAllFreqMean, [2 1 3])); % NOTE: spectral data is exponentiated in gamma_spectral_analysis

figure; set(gcf, 'color', 'w');
for chan = 1:size(specData,3)
    subplot(1,2,1); cla
    set(gca, 'Colororder', colormap, 'XScale', 'log', 'XLim', [35 200]); hold on
    plot(f_plot, specData(:, :, chan),...
        'LineWidth', 2);
    legend(conditionNames)
    yl = get(gca, 'YLim');
    
    subplot(1,2,2); cla
    set(gca, 'Colororder', colormap, 'XScale', 'log', 'XLim', [35 200]); hold on
    plot(f_plot, fit(:, :, chan),...
        'LineWidth', 2);
    set(gca, 'YLim', yl);
    pause(1)
    if opt.saveFigures; if ~exist(fullfile(opt.sessionPath, 'figs', 'spectraCondChan'),'dir'); 
            mkdir(fullfile(opt.sessionPath, 'figs', 'spectraCondChan')); end
            hgexport(gcf, fullfile(opt.sessionPath, 'figs', 'spectraCondChan',saveName)); end
end

%% 2. All conditions' line fit plotted on the same spectrogram for each channel

figure; set(gcf, 'color', 'w');
for chan = 1:size(specData,3)
    clf; set(gcf, 'Name', sprintf('Channel %d', chan));
    hold all;
    
    for ii = 1:length(conditionNames)
        plot(f_plot, smooth(fit(:, ii, chan), 2)',...
            'color', colormap(ii,:,:),'LineWidth', 2);
    end
    legend(conditionNames)
    title(sprintf('Modelfit for channel %d', chan))
    pause(1);
    if opt.saveFigures; if ~exist(fullfile(opt.sessionPath, 'figs', 'spectraSmoothed'),'dir'); 
            mkdir(fullfile(opt.sessionPath, 'figs', 'spectraSmoothed')); end
            hgexport(gcf, fullfile(opt.sessionPath, 'figs', 'spectraSmoothed',saveName)); end
end

%% All channels, fitted spectrogram for each condition

clf;
figure(3);
for chan = 1:size(specData,3)
    set(gcf, 'Name', sprintf('Channel %d', chan));
   for cond = 1:length(conditionNames)-1
       subplot(4,3,cond); cla;
       % plot given stimuli condition
       plot(f_plot, fit(:,cond,chan), 'Color', colormap(cond,:)); hold on;
       % plot baseline fit
       plot(f_plot, fit(:,length(conditionNames), chan), 'Color', 'k');
      
       %set(gca, 'Color', colormap(cond, :));
       title(cell2mat(conditionNames(cond)));
       title(sprintf('%s BB:%f1.3;  G:%f1.3', cell2mat(conditionNames(cond)), results.broadbandPowerMean(chan, cond), results.gammaPowerMean(chan, cond)));
       set(gca, 'XScale', 'log', 'XLim', [35 200]);
   end
   pause(1)
   if opt.saveFigures; if ~exist(fullfile(opt.sessionPath, 'figs', 'spectraSeparateConds'),'dir');
            mkdir(fullfile(opt.sessionPath, 'figs', 'spectraSeparateConds')); end
            hgexport(gcf, fullfile(opt.sessionPath, 'figs', 'spectraSeparateConds',saveName)); end

end

