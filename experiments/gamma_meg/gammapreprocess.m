function newdata = gammapreprocess(data, t, f, evoked_cutoff, keep_frequencies)
% Clip and filter time series prior to denoising for gamma experiment. We
% clip the early part of the epoch to remove the evoked response. From the
% remaining response, we filter out all frequencies that are not used for
% extracting the broadband and gamma responses.
%
% data = gammapreprocess(data, t, f, keep_timepts, keep_frequencies)

% data is [channels x time points x epochs]

% clip the data and time vector
data    = data(:, t>=evoked_cutoff, :);
ln      = 60; % line noise

% harmonics of of the stimulus locked
tmp = (1:10) * ln; 
drop_frequencies  = [f(sort(unique([tmp-1 tmp tmp+1]))), ~keep_frequencies];
% high pass filter with cutoff of 62 Hz, sharp cutoff, and excluding
% harmonics
fs = 1000;
newdata = filterdata(data,fs,0,1,drop_frequencies);


return