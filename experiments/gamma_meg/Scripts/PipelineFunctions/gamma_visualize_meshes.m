function [] = gamma_visualize_meshes(results, opt)
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

%% Add fieldtrip path if necessary
if isempty(which('ft_prepare_layout'))
    addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'))
end

% In case data is analyzed on the HPC
opt.sessionPath = meg_gamma_get_path(opt.params.sessionNumber);
opt.saveFigures = 1;

%% Prepare data
conditionNames = gamma_get_condition_names(opt.params.sessionNumber);

% Define the file name of saved figures
if opt.saveFigures
    if opt.MEGDenoise; postFix = '_denoised'; else postFix = ''; end
    preFix = sprintf('s_%3d', opt.params.sessionNumber);
    saveName = sprintf('%s_boots%d%s_%s', preFix, opt.nBoot, postFix,...
        datestr(now, 'mm_dd_yy'));
end

%% Define summary metric (SNR or just mean)

if opt.nBoot > 1,
    summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
else
    summary_stat = @(x) nanmean(x,3);
end

numConds = length(unique(results.opt.params.conditions));
numChannels = size(results.spectralData,3);

% Get difference identity vectors (each condition - baseline)
contrasts = eye(numConds,numConds);
% Make baseline column negative to subtract from each cond
contrasts(:,results.opt.params.baselineCondition) = -1;

contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
numContrasts = size(contrasts,1);

% compute SNR
tmpData = permute(results.broadbandPower, [2 1 3]);
tmpData = reshape(tmpData, numConds, []);
tmp = contrasts*tmpData;
tmp = reshape(tmp, numContrasts, numChannels, opt.nBoot);
snrBroadband = summary_stat(tmp)';

tmpData = permute(results.gammaPower, [2 1 3]);
tmpData = reshape(tmpData, numConds, []);
tmp = contrasts*tmpData;
tmp = reshape(tmp, numContrasts, numChannels, opt.nBoot);
snrGamma  = summary_stat(tmp)';

%% Plot meshes Gamma

scrsz = get(0,'ScreenSize');
threshold = 0;%3;

figure; clf; set(gcf, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)]);  set(gcf, 'name', 'Gaussian SNR' )
plotRange = [-1 1] * (max(max(abs(snrGamma(:,1:numConds)))));
for c = 1:12
    subplot(4,3,c)
    dataToPlot = snrGamma(:,c)';
    dataToPlot(abs(dataToPlot) < threshold) = 0;
    if length(dataToPlot) < 157; ft_plotOnMesh(to157chan(dataToPlot,opt.params.badChannels,'nans'), conditionNames{c});
    else ft_plotOnMesh(dataToPlot, conditionNames{c}); end
    set(gca, 'CLim', plotRange)
    colormap bipolar
end

if opt.saveFigures;
    if ~exist(fullfile(opt.sessionPath,'figs','meshes'),'dir'); mkdir(fullfile(opt.sessionPath,'figs','meshes')); end
    hgexport(gcf, fullfile(opt.sessionPath,'figs','meshes',sprintf('SNR_Gamma_mesh_%s_%d.eps', saveName, threshold)));
end


%%
%SNR Power

scrsz = get(0,'ScreenSize');
threshold = 0;%3;
figure; clf, set(gcf, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)]); set(gcf, 'name', 'Broadband SNR')
plotRange = [-1 1] * (max(max(abs(snrBroadband(:,1:numConds)))));

for c = 1:12
    subplot(4,3,c)
    dataToPlot = snrBroadband(:,c)';
    dataToPlot(abs(dataToPlot) < threshold) = 0;
    if length(dataToPlot) < 157; ft_plotOnMesh(to157chan(dataToPlot,opt.params.badChannels,'nans'), conditionNames{c});
    else ft_plotOnMesh(dataToPlot, conditionNames{c}); end
    set(gca, 'CLim', plotRange)
    colormap bipolar
end

if opt.saveFigures;
    if ~exist(fullfile(opt.sessionPath,'figs','meshes'),'dir'); mkdir(fullfile(opt.sessionPath,'figs','meshes')); end
    hgexport(gcf, fullfile(opt.sessionPath,'figs','meshes',sprintf('SNR_Broadband_mesh_%s_%d.eps', saveName, threshold)));
end


