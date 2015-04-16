function epoch_ts = make_epoch_ts(conditions, onsets, n_samples)
%
% Create a timeseries that marks the start of every epoch, as opposed to
% every trigger. Additionally, meg_make_epochs requires a ts with
% stimulus condition numbers (e.g. 1, 3, 5, and 7 for full visual field,
% blank, left field, right field in the SSEEG experiment), instead of just
% using ones to mark events, as we did in ev_ts. 
%   
% INPUTS: 
%    condtions:     A cell (1 x nr_runs) of vectors (
%    onsets:        A cell of vectors representing the frame number that each
%                   epochs starts
%    n_samples:     A vector indicating the number of samples in each run.
%                   If a scalar, then assume the same number for every run
%
% OUTPUTS
%   epoch_ts:       a 
nr_runs = numel(conditions);

% one condition vector and one onsets vector for each run. so they
% must be the same length
try
    assert(numel(conditions) == numel(onsets))
catch
    error('The number of runs with conditions is not the same as the number of runs with onset times')
end


if isscalar(n_samples), n_samples = repmat(n_samples, [1, nr_runs]); end

epoch_ts = cell(1, nr_runs);

for ii = 1:nr_runs
    epoch_ts{ii} = zeros(1,n_samples(ii));    
    epoch_ts{ii}(onsets{ii}) = conditions{ii};    
end

return
