function [sensorData, badChannels, badEpochs] = meg_preprocess_data(sensorDataIn, varThreshold, ...
    badChannelThreshold, badEpochThreshold, sensorPositions, verbose)
% Preprocess MEG data
%
%sensorData = meg_preprocess_data(sensorDataIn, varThreshold, ...
%    badChannelThreshold, badEpochThreshold, sensorPositions, verbose)
%
% INPUTS
%   sensorDataIn: 3D array, time points x epochs x channels
%
%   varThreshold: Vector of length 2 ([min max]) to indicate variance
%                   threshold. For any channel in any give epoch, if the
%                   variance in the time series is outside [min max] * the
%                   median variance across all channels and all epochs,
%                   then label that channel during that epoch as 'bad'
%                       Default = [.05 20]
%
%  badChannelThreshold: Fraction ([0 1]). If more than this fraction of
%                   channels in any given epoch is labeled 'bad', then
%                   label all channels for this epoch as 'bad'.
%                       Default = 0.2
%
%  badEpochThreshold: Fraction ([0 1]). If more than this fraction of
%                   epochs for any given channel is labeled 'bad', then
%                   label all epochs for this  channel as 'bad'
%                       Default = 0.2
%  sensorPositions: mat file name with xyz channel positions OR matrix with
%                   xyz positions (n channels x 3)
%
%  verbose:          Whether to plot debug figures and display info
%
% OUTPUTS
%   sensorData:     Same as sensorDataIn (3D array, time points x epochs x
%                   channels), except that for which 'bad', data has been
%                   replaced with NaN
%
%
% Example:
%  sensorData = dfdPreprocessData(sensorDataIn, [.01 10], .2, .2);

if notDefined('varThreshold'), varThreshold = [.05 20]; end
if notDefined('badChannelThreshold'), badChannelThreshold = .2; end
if notDefined('badEpochThreshold'), badEpochThreshold = .2; end
if notDefined('sensorPositions'), sensorPositions = 'meg160xyz'; end
if notDefined('verbose'), verbose = false; end

% This identifies any epochs whos variance is outside some multiple of the
% grand variance
outliers = meg_find_bad_epochs(sensorDataIn, varThreshold);

% any epoch in which more than 10% of channels were bad should be removed
% entirely
badEpochs = mean(outliers,2)>badChannelThreshold;

% once we remove 'badEpochs', check whether any channels have more
% than 10% bad epochs, and we will remove these
badChannels = mean(outliers(~badEpochs,:),1)>badEpochThreshold;

outliers(badEpochs,:)   = 1;
outliers(:,badChannels) = 1;

% Plot outiers for epochs and channels
if verbose
    figure; imagesc(outliers);
    xlabel('channel number'); ylabel('epoch number'); title('Bad channels / epochs')
    fprintf('(dfdPreprocessData): %5.2f%% of epochs removed\n', sum(sum(outliers))/(size(sensorDataIn,2)*size(sensorDataIn,3))*100);
end

% Interpolate epochs over neighbouring channels
sensorData = meg_channel_repair(sensorDataIn, outliers, 'nearest', sensorPositions);


return


