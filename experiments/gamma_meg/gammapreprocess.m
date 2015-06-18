function data = gammapreprocess(data, t, f, evoked_cutoff, keep_frequencies)
% Clip and filter time series prior to denoising for gamma experiment. We
% clip the early part of the epoch to remove the evoked response. From the
% remaining response, we filter out all frequencies that are not used for
% extracting the broadband and gamma responses.
%
% data = gammapreprocess(data, t, f, keep_timepts, keep_frequencies)

% data is [channels x time points x epochs]

% clip the data and time vector
data    = data(:, t>evoked_cutoff, :);
t       = t(t>evoked_cutoff);

error('Not yet finished')