%% s_GAMMA_ChannelSelection

% We want a function that shows which channels which show highest activity
% based on different criteria (type of signal, image, clusters?)

% Workflow
% 0. Get session information
% 1. Define summary metric for modelfits used for channel selection
% 2. Select channels
% 3. Plot mesh showing the selected channels
% 4. Plot the summarized spectra

%% 0. Get session information
sessionNum      = 20;
sessionPath     = meg_gamma_get_path(sessionNum);
conditionNames  = gamma_get_condition_names(sessionNum);

% Get sessions parameters
params          = gamma_get_parameters(sessionNum);
opt.params      = params; clear params;
opt.sessionPath = sessionPath; clear sessionPath;
opt.saveFigures = 0;

% Selection method
selectionMethod = 'gamma'; % use gamma
summaryStat     = @(x) nanmean(x,3) ./ nanstd(x, [], 3); % use SNR
topX            = 10;

% Define the file name of saved figures
if opt.saveFigures
    if opt.MEGDenoise; postFix = '_denoised'; else postFix = ''; end
    preFix = sprintf('s_%3d', opt.params.sessionNumber);
    saveName = sprintf('%s_boots%d%s_%s', preFix, opt.params.nBoot, postFix,...
        datestr(now, 'mm_dd_yy'));
end

% Load data
d = dir(fullfile(opt.sessionPath,'processed',sprintf('s0%02d_results_100boots_denoised_*',sessionNum)));
tmp = load(fullfile(opt.sessionPath,'processed',d.name)); results = tmp.results;


%%  1. Define summary metric (SNR or just mean of modelfits)

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
tmp = reshape(tmp, numContrasts, numChannels, opt.params.nBoot);
snrBroadband = summaryStat(tmp)';

tmpData = permute(results.gammaPower, [2 1 3]);
tmpData = reshape(tmpData, numConds, []);
tmp = contrasts*tmpData;
tmp = reshape(tmp, numContrasts, numChannels, opt.params.nBoot);
snrGamma  = summaryStat(tmp)';


%%  2. Select top 5 and 10 SNR channels from stimulus conditions for the two
% types of signals



switch selectionMethod
    case 'gamma'
        nanChan = find(isnan(snrGamma(:,1)));
        snrGamma(nanChan,:) = -inf;
        snrGamma(results.opt.params.badChannels,:) = -inf;
        [value, idx] = sort(max(snrGamma,[],2),'ascend');
end

%%  3. Show top X selected channels on a Mesh (not implemented yet)
% cfg.highlight          = 'on'; % or use 'labels', 'numbers', 'off'
% cfg.highlightchannel   = idx(1:topX)'; % Nx1 cell-array with selection of channels, or vector containing channel indices see FT_CHANNELSELECTION
% cfg.highlightsymbol    = 'o'; %highlight marker symbol (default = 'o')
% cfg.highlightcolor     = [0 0 0];%highlight marker color (default = [0 0 0] (black))
% cfg.highlightsize      = 12; %highlight marker size (default = 6)
% cfg.highlightfontsize  = 12; %highlight marker size (default = 8)

% ft_plotOnMesh(snrGamma(:,1),conditionNames{1}, [], '2d', cfg)

%% Plot spectra

% Get spectra
% d = dir(fullfile(opt.sessionPath,'processed',sprintf('s0%02d_spectralData_100boots_denoised_*',sessionNum)));
d = dir(fullfile(opt.sessionPath,'processed',sprintf('s0%02d_spectralData_100boots_07*',sessionNum)));

tmp = load(fullfile(opt.sessionPath,'processed',d.name)); spectralDataBoots = tmp.spectralDataBoots;

spectralData = mean(spectralDataBoots,4);

%% If you want to manually change opt.fitFreq --
numFreq        = size(spectralData, 1); % the number of frequency bins in data
t              = (1:numFreq)/results.opt.fs;
f              = (0:length(t)-1)/max(t);

% fitFreq = f((f>=35 & f <= 59) | (f>=61 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
fitFreq = f((f>=35 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));


f_sel          = ismember(f, fitFreq);
f_plot         = f;
f_plot(~f_sel) = NaN;

% spectra sorted by condition..
figure;
for ii = 1:12
    subplot(4,3,ii);
    plot(f_plot, nanmean(squeeze(spectralData(:,ii,idx(1:topX))),2), 'LineWidth',4); hold all;
    plot(f_plot, nanmean(squeeze(spectralData(:,results.opt.params.baselineCondition,idx(1:topX))),2), 'k')

        
    for jj = 1:topX
        plot(f_plot, nanmean(squeeze(spectralData(:,ii,idx(jj))),2)); end
    
    title(conditionNames(ii))
    
    set(gca, 'XScale', 'log', 'XLim', [35 200]);
    xlabel('Frequency (Hz)')
    ylabel('Power (fT.^2)')
end











