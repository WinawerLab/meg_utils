function [meanResults, opt] = gamma_spectral_analysis(ts, opt)
%% results = gamma_spectral_analysis(ts, opt)
%  performs the spectral analysis, bootstrapping, and curve fitting for
%  s_GAMMA_pipeline

% INPUTS
% ts   : epoched timeseries, can be denoised or not yet denoised (time x epochs x channels)
%      : opt - structure containing experimental parameters obtained
% 
% OUTPUTS
% results (a struct which includes):
%       spectralDataMean    :    3D array containing spectral data averaged across bootstraps (time x conditions x channels) 
%       fitBaselineMean     :    2D array containing baseline fit averaged across bootstraps (channels x selected frequencies) 
%       broadbandPowerMean  :    2D array containing broadband power elevation fit averaged across bootstraps (channels x conditions) 
%       gammaPowerMean      :    2D array containing narrowband gamma power gaussian fit averaged across bootstraps (channels x conditions) 
%       gammaPeakFreqMean   :    1D array containing peak frequency of gaussian bump (channels x 1) 
%       modelFitAllFreqMean :    3D array containing modelfit for each frequency, condition, channel (frequencies, conditions, channels)
%       fitFreq             :    1D array containing all frequencies used for model fit
%       opt                 :    struct with options used when analyzing this dataset

%% Get postFix for this spectral analysis based on parameters
thisDate = datestr(now, 'mm.dd.yy');
if opt.MEGDenoise; denoiseStr = '_denoised'; else denoiseStr = ''; end
postFix   = sprintf('%dboots%s_%s', opt.nBoot, denoiseStr, thisDate);

%% Calculate Spectral Data

% Get conditions
conditionsUnique = unique(opt.params.conditions);
numConditions = length(conditionsUnique);

% Parameters for spectral analysis and modelfitting
t = (1:size(ts,1))/opt.fs;
f = (0:length(t)-1)/max(t);
opt.fitFreq = f((f>=35 & f <= 56) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200)); % To do: get this hard coded part out of the function 

% FFT over time points
spectralData = abs(fft(ts))/length(t)*2;

% Prepare array for bootstrapping
% (freq x conditions x channels x bootstraps)
spectralDataBoots = zeros(size(ts,1), length(conditionsUnique), length(opt.dataChannels), opt.nBoot);

% compute the mean amplitude spectrum for each electrode in each condition
if opt.verbose; fprintf('[%s]: Computing bootstraps for each condition\n', mfilename); end
for ii = 1:length(conditionsUnique)
    if opt.verbose; fprintf('[%s]: Condition %d of %d\n', mfilename, ii, length(conditionsUnique)); drawnow; end
    
    % Binary vector to identify epochs with this condition
    theseEpochs = opt.params.conditions == conditionsUnique(ii);
    
    % Get spectral data of one condition (time points x epochs x channel)
    theseData = spectralData(:,theseEpochs,opt.dataChannels);
    
    if opt.nBoot > 1
        % Reshape data so that epochs are in rows (needed for bootstrapping)
        theseData = permute(theseData, [2 1 3]);
        
        % Take the log normalized mean of spectral power
        bootfun = @(x) squeeze(exp(nanmean(log(x),1)));
        
        % Bootstrap (this gives us a matrix: nboot x (freq x channel)
        bootstat = bootstrp(opt.nBoot, bootfun, theseData);
        
        % Reshape bootstatistics into 3D-array: nboot x freq x channel
        bootstat = reshape(bootstat, opt.nBoot, length(t), []);
        
        % Now permute the data so we get spectralDataBoots is freq x condition x channel x boot
        spectralDataBoots(:,ii,:,:) = permute(bootstat,[2 3 1]);   
        
    else % If no bootstrapping
        spectralDataBoots(:,ii,:,:) = exp(nanmean(log(theseData),2));
    end
end
if opt.verbose; fprintf('[%s]: Done!\n', mfilename); end

% Summarize bootstrapped spectral by mean and std over bootstraps
spectralDataMean = nanmean(spectralDataBoots, 4);


%% Fit broadband and Gaussian power by a model that exists of an parallel elevation from baseline and a gaussian bump

% Binary array describing the frequiencies not ommitted from the fit
opt.params.f_sel = ismember(f,opt.fitFreq);

% per loop sizes
% fit_bl: chan x 112 x boot
% w_pwr: chan x cond x boot
% gauss_f : chan x 1 x boot
% fit_f2: freq x cond x chan x boot

numChannels      = size(spectralDataMean,3);
fitBaseline      = NaN(numChannels,length(opt.fitFreq), opt.nBoot);   % slope of spectrum in log/log space
broadbandPower   = NaN(numChannels,numConditions, opt.nBoot);     % broadband power
gammaPower       = NaN(numChannels,numConditions, opt.nBoot);     % gaussian height
gammaPeakFreq    = NaN(numChannels, 1, opt.nBoot);                 % gaussian peak frequency
modelFitAllFreq  = NaN(numConditions,size(spectralDataMean,1),numChannels, opt.nBoot); % fitted spectrum

warning off 'MATLAB:subsassigndimmismatch'

% For each channel, fit each condition separatley
fprintf('[%s]: Fitting gamma and broadband values for each channel and each condition', mfilename)

% Fit each channel separately
for chan = 1:size(spectralDataMean,3)
    
    if opt.verbose; fprintf('[%s]: Channel %d of %d\n', mfilename, chan, size(spectralDataMean,3)); drawnow; end
    
    % Fit each bootstrap separately   
    for bootnum = 1:opt.nBoot
        
        % the baseline is the same for all conditions, so just compute it once
        dataBase = exp(mean(log(spectralDataBoots(:,:,chan, bootnum)),2));
        dataBase = dataBase';
        
        if all(isfinite(dataBase)) && mean(dataBase) > 0           
            dataFit = spectralDataBoots(:,:,chan, bootnum);            
            [...
                fitBaseline(chan, :, bootnum), ...
                broadbandPower(chan, :, bootnum), ...
                gammaPower(chan, :, bootnum),...
                gammaPeakFreq(chan, 1, bootnum),...
                modelFitAllFreq(:,:, chan, bootnum)] = ...
                gamma_fit_data_localregression_multi(f,opt.fitFreq,dataBase,dataFit);
        end
    end
    
    % For debugging:
    if opt.verbose
    figure(1); clf;
    plot(f, spectralDataBoots(:,:,chan, bootnum), 'k'); hold all;
    plot(f, exp(modelFitAllFreq(:,:, chan, bootnum')),'r')
    set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [10 200], 'YLim', 10.^[0 3])
    title(sprintf('Channel %d', chan))
    pause(0.1)
    if opt.saveFigures; hgexport(fullfile(meg_gamma_get_path(sessionNum),'figs',sprintf('spectrumAllCondsChan%d',chan))); end
    end
end

if opt.verbose; fprintf('[%s]: done!\n', mfilename); end

warning on 'MATLAB:subsassigndimmismatch'

% Summarize bootstrapped fits
gammaPeakFreq        = repmat(gammaPeakFreq, [1 numConditions 1]);
fitBaselineMean      = nanmean(fitBaseline,3);
broadbandPowerMean   = nanmean(broadbandPower,3);
gammaPowerMean       = nanmean(gammaPower,3);
gammaPeakFreqMean    = nanmean(gammaPeakFreq,2);
modelFitAllFreqMean  = nanmean(modelFitAllFreq,4);

results = struct('spectralData', spectralDataMean, ...
                 'fitBaseline', fitBaselineMean, ...
                 'broadbandPower', broadbandPowerMean, ...
                 'gammaPower', gammaPowerMean, ...
                 'gammaPeakFreq', gammaPeakFreqMean, ...
                 'modelFitAllFreq', modelFitAllFreqMean, ...
                 'fitFreq', opt.fitFreq, ...
                 'opt', opt);

% Summarize all in one results struct
meanResults = struct('spectralDataMean', spectralDataMean, ...
                 'fitBaselineMean', fitBaselineMean, ...
                 'broadbandPowerMean', broadbandPowerMean, ...
                 'gammaPowerMean', gammaPowerMean, ...
                 'gammaPeakFreqMean', gammaPeakFreqMean, ...
                 'modelFitAllFreqMean', modelFitAllFreqMean, ...
                 'fitFreq', opt.fitFreq, ...
                 'opt', opt);
             


%% Save spectral data and results

if opt.saveData 
    if opt.verbose; fprintf('[%s]: Saving data..]',mfilename); end
    if opt.HPC
        savePth = fullfile(opt.sessionPath, 'processed');
    else
    savePth = fullfile(meg_gamma_get_path(opt.params.sessionNumber), 'processed'); end
    % Save bootstrapped spectral data
    save(fullfile(savePth,sprintf('s%03d_spectralData_%s.mat', opt.params.sessionNumber, ...
        postFix)),'spectralDataBoots', '-v7.3')
    % Save results
    save(fullfile(savePth,sprintf('s%03d_results_%s.mat', opt.params.sessionNumber, ...
        postFix)),'results', '-v7.3')
    % Save mean results
    save(fullfile(savePth,sprintf('s%03d_meanresults_%s.mat', opt.params.sessionNumber, ...
        postFix)),'meanResults', '-v7.3')
end

end