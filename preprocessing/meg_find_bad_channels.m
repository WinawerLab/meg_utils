function bad_channels = meg_find_bad_channels(ts)
%Identify problem channels in MEG time series. 
%
%  bad_channels = meg_find_bad_channels(ts)
%
%   Bad channels are identified based on the standard deviation of the
%   time series. We compute one sd for each channel for each epoch. Then,
%   for each channel, we take the median standard deviation across epochs.
%   The grand median is the median of these values across channels. If the
%   median of any particular channel is either less than 10 times smaller
%   or more than 10 times bigger than the grand median, we label that
%   channel as bad. We return a vector of bad channels.
% 
% INPUT
%   ts: time x epoch x channel
% OUTPUT
%   bad_channels: vector of bad channels
%
% Example: bad_channels = meg_find_bad_channels(ts);



sd              = squeeze(nanstd(ts)); % sd will be epochs x channels
sd_median       = nanmedian(sd); % median sd for each channel
sd_grand_median = nanmedian(sd_median); % median across all channels

bad_channels = sd_median < 0.1 * sd_grand_median | sd_median > 10 * sd_grand_median;
bad_channels = find(bad_channels);




