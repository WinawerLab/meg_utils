function evoked = getevoked(data,fs, design_mtrx, peak_window)
% get the amplitude of the evoked signal
%
% evoked = getevoked(data,opts)
%
% INPUT
%    data:          raw time series [channels x time x epochs]
%
%    fs:            sampling frequency (in Hz) [default = 1000]
%
%    design_mtrx:   binary matrix, epochs x conditions [dafault = ones(num_epochs, 1)
%
%    peak_window:   [start end] relative to grand peak (in samples):
%                       This determines the possible time points in which
%                       we will accept a peak for individual channels.
%                       [default = [-30 30]
%
% OUTPUT
%    evoked : evoked time series, 1 number for each of [epochs x channels]


% check input 
if notDefined('fs'),     fs     = 1000; end
if notDefined('design_mtrx'), design_mtrx = ones(size(data, 2), 1); end
if notDefined('window'), peak_window = [-30 30]; end

% reshape data
data         = permute(data, [2 3 1]); % [time x epochs x channels]
num_epochs   = size(data, 2);
num_channels = size(data, 3);

% check that size of data and size of design matrix match
assert(num_epochs == size(design_mtrx, 1));

% low pass filter the data
df = designfilt('lowpassfir','FilterOrder',70,'CutoffFrequency',50, 'SampleRate', fs);
D  = mean(grpdelay(df)); % filter delay in samples
ts_f = padarray(data, [D 0 0], 0, 'post'); % Append D zeros to the input data
ts_f = filter(df,ts_f);                    % Filter it!
ts_f = ts_f(D+1:end,:, :);                 % Shift data to compensate for delay

% don't include blank epochs in finding grand peak
blank_trials = sum(design_mtrx, 2) == 0;

% To find the peak evoked time point, we do two steps:
%       1. Compute the average temporal response of each channel across all
%           non-blank epochs
evoked_grand     = squeeze(nanmean(ts_f(:, ~blank_trials, :) , 2));
%       2. Compute the rms at each time point across channels
evoked_grand_rms = sqrt(nansum(evoked_grand.^2,2));
% grand_idx is the time point (in samples) with the highest RMS
[~, grand_idx] = max(evoked_grand_rms);

% We expand the seach window so that the peak of each channel can vary a bit 
peak_window = peak_window + grand_idx;

% We clip the window if it is bigger than the actual epoch length
if peak_window(2) > size(data,1)
    peak_window(2) = size(data,1);
end

if peak_window(1) < 1, peak_window(1) = 1; end

% which_time_point is a vector of length num_channels that tells us which
% time point is the peak for each channel
% [~, which_time_point] = max(abs(evoked_grand(peak_window, :)));
[~, which_time_point] = max(abs(evoked_grand(min(peak_window):max(peak_window), :)));
which_time_point = which_time_point + min(peak_window) -1;

% Finally, we go back to the data, and pull out the value at the
% appropriate time point for each channel from each epoch and thereby
% return a matrix with one number per channel per epoch
evoked = NaN(num_epochs, num_channels);

for ii = 1:num_channels
  evoked(:,ii) = data(which_time_point(ii), :, ii);
end

return