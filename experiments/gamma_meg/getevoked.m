function evoked = getevoked(data,opts)
% get the amplitude of the evoked signal
% INPUT
% data   : raw time series [channels x time x epochs]
% opts   : 
%           fs: sampling frequency (in Hz)
%           design: binary matrix, epochs x conditions
% OUTPUT
% evoked     : evoked time series, 1 number for each of [epochs x channels]

% check input 

data = permute(data, [2 3 1]); % [time x epochs x channels]

num_epochs = size(data, 2);
num_channels = size(data, 3);

df = designfilt('lowpassfir','FilterOrder',70,'CutoffFrequency',50, 'SampleRate', opts.fs);
D  = mean(grpdelay(df)); % filter delay in samples

ts_f = padarray(data, [D 0 0], 0, 'post'); % Append D zeros to the input data
ts_f = filter(df,ts_f);                    % Filter it!
ts_f = ts_f(D+1:end,:, :);                 % Shift data to compensate for delay

blank_trials = sum(opts.design, 2) == 0;

evoked_grand     = squeeze(nanmean(ts_f(:, ~blank_trials, :) , 2));

evoked_grand_rms = sqrt(nansum(evoked_grand.^2,2));
[~, grand_idx] = max(evoked_grand_rms);
peak_window = (-30:30) + grand_idx;

[~, which_time_point] = max(abs(evoked_grand(peak_window, :)));
which_time_point = which_time_point + min(peak_window) -1;

evoked = NaN(num_epochs, num_channels);

for ii = 1:num_channels
  evoked(:,ii) = data(which_time_point(ii), :, ii);
end