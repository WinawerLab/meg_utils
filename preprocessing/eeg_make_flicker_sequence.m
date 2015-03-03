function flicker_sequence = eeg_make_flicker_sequence(n, dur, isi, s_rate)
error('Not yet implemented')
%
% make a sequence of n flashes which each last for dur seconds, separated
% by isi seconds, sampled at s_rate
%
% Example: 4 50-ms flases with 100 ms spacing at 1000 Hz recording
% n = 4;
% dur = .050
% isi = .100
% s_rate = 1000;
% flicker_sequence = eeg_make_flicker_sequence(n, dur, isi, s_rate);

flicker_sequence = zeros(1, round(n*(dur+isi)/dt));

