function ts_denoised = meg_environmental_denoising(ts, opt)
% This function uses the timeseries of the three MEG non-physiological
% channels to denoise the physiological channels. 
%  
% ts_denoised = meg_environmental_denoising(ts,...
%     environmental_channels, data_channels, ...
%     [produce_figures=0], [save_data=0], [verbose=0])
%
% INPUTS:
%  ts                        MEG timeseries (number of timepoints by number 
%                                of epochs  by number of channels)
%  environmental_channels    vector of channel numbers that are
%                               nonphysiological. These will be used to
%                               regress out noise from all other channels.
%  data_channels             vector of channel numbers to denoise
%  produce_figures           boolean. If true, make some plots to compare
%                               pre and post denoising
%  verbose
%
% OUTPUTS:
%  ts_denoised               denoised time series


%% Make empty arrays for regressed 'clean' data
ts_denoised = ts;

% Start regression, keep residuals
for epoch = 1:size(ts,2); % Epoch size is the same for every condition (i.e. 180 except for session 3 (=168)) <--- EK: so do we need to do anything with this?!
    projectOut  = squeeze(ts(:,epoch,opt.environmentalChannels));
    projectFrom = squeeze(ts(:,epoch,opt.dataChannels));
    projectionWeights  = projectOut \ projectFrom;
    ts_denoised(:,epoch,opt.dataChannels) = projectFrom - projectOut*projectionWeights;
end

%% For debugging: Make figures of all the raw epochs 
if opt.verbose
    % And for two visual channels (channel 1 and 14)
    figure; plot(squeeze(nanmean(ts_denoised(:,:,1),2)),'r'); hold on;
    plot(squeeze(nanmean(ts(:,:,1),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (pTesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', 1))
    legend('Denoised','Raw')

    figure; plot(squeeze(nanmean(ts_denoised(:,:,14),2)),'r'); hold on;
    plot(squeeze(nanmean(ts(:,:,14),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (pTesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', 14))
    legend('Denoised','Raw')
end


