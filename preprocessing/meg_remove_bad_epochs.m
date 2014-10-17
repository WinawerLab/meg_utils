function ts = meg_remove_bad_epochs(bad_epochs, ts)

% This function will take every epoch and compares the variance of the
% signal within this epoch with a certain threshold (say, 10 times the
% median of the signal of one channel). All epochs with a variance larger
% than the threshold, will be defined as a bad epoch (1) in a matrix. One 
% can use this matrix to plot the difference between with/without bad
% epochs.

% INPUTS
%   bad_epochs:       
%   ts:           Timeseries (time points x epochs x channels)
%
% OUTPUTS 
%   ts:         Timeseries with removed epochs


%% Make matrix to remove bad epochs for ON-data

% epochs x channel
num_time_points = size(ts,1);
num_epochs      = size(ts, 2);
num_channels    = size(ts, 3);

ts = reshape(ts, [num_time_points, num_epochs*num_channels]);

ts(:, logical(bad_epochs(:))) = NaN;

ts = reshape(ts, [num_time_points, num_epochs, num_channels]);

return
