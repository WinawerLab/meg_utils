function ts = meg_remove_bad_channels(bad_channels, ts)
%Replace data for bad channels with NaN
%
% ts  = meg_remove_bad_channels(bad_channels, ts)
%
% INPUTS:
%   bad_channels: vector of bad channels
%   ts: data in format time x epochs x channels
%
% OUTUTS
%   ts: same as input ts but with data from bad channels replaced with NaN
%
% Example: ts = ssmeg_remove_bad_channels([137 153], ts);

ts(:,:, bad_channels) = NaN;

return
