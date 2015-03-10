function flicker_seq = eeg_make_flicker_sequence(n, dur, isi, s_rate, pad_with_zeros)
%
% flicker_sequence = eeg_make_flicker_sequence(n, dur, isi, s_rate, pad_with_zeros);
% 
% Function to make a sequence of 1's and 0's to define the order of the 
% photodiode flashes with duration, and isi in seconds, sampled at s_rate. 
% 
% Inputs: 
% n     :   Number of flashes
% dur   :   Duration of a flash in seconds
% isi   :   Inter-stimulus-interval, time in seconds in between flashes
% s_rate:   Sample rate of photo-diode data, in Hz
%
% Example: 4 50-ms flashes with 100 ms spacing at 1000 Hz recording
% n = 4;
% dur = .050
% isi = .100
% s_rate = 1000;
% 
% Get the length of the sequence in time points
n_timepoints = (n*dur + (n-1)*isi) * s_rate;

% Make an empty flicker sequence
flicker_seq = zeros(1,round(n_timepoints)-1);

% Fill it up
for ii = 1:n
    flicker_seq(round(((ii-1)*dur*s_rate) + ((ii-1)*isi*s_rate) + (1:(dur*s_rate)) ))=1;
end

% pad with zeros
if exist('pad_with_zeros', 'var')
    flicker_seq = padarray(flicker_seq(:), pad_with_zeros);
end

return