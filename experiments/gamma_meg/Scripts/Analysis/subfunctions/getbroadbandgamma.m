function ab = getbroadbandgamma(data,freq)
% get the amplitude of the evoked signal
%
% ab = broadband_or_gamma(data,freq)
%
% INPUT
%    data   : raw time series [channels x time x epochs]
%    opts   :  fs:     sampling frequency (in Hz) [default = 1000]
%              design: binary matrix, epochs x conditions [dafault = ones(num_epochs, 1)
%              window: [start end] relative to grand peak (in samples):
%                       This determines the possible time points in which
%                       we will accept a peak for individual channels.
%                       [default = [-30 30]
% OUTPUT
%    evoked : broadband_or_gamma time series, 1 number for each of [epochs x channels]


% check input 
if isnumeric(freq)
    f = freq;
elseif isstruct(freq) && isfield(freq,'ab_i')
    f = freq.ab_i;
else
    error('input error: freq not recognized');
end

spec = fft(data,[],2);
spec_amp = abs(spec(:,f,:))/ size(data,2)*2; 
% take the mean across frequencies in log space 
ab = squeeze(exp(nanmean(log(spec_amp.^2),2)))';

