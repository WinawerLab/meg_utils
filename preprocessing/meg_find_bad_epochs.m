function bad_epochs = meg_find_bad_epochs(ts, thresh)
%Identify problem epochs in MEG time series. 
%
%  bad_epochs = meg_find_bad_epochs(ts)
%
%   Bad epochs are identified based on the variance of the
%   time series. We compute the variance for each channel for each
%   epoch. Any of these that are above some multiple (default = 10) of the
%   grand median, or below some fraction (default = 1/10) of the grand
%   median are considered outliers.
% 
% INPUT
%   ts:     time x epoch x channel
%   thresh: 2-vector specificy min and max multiplier of the variance to
%                 threshold data
% OUTPUT
%   bad_epochs: logical matrix of bad epocjs (epoch x channel), where 1 =
%                   bad
%
% Example: bad_epochs = meg_find_bad_channels(ts, [.1 10]);

if ~exist('thresh', 'var') || isempty(thresh), 
    thresh = [0.1 10];
end

var_matrix       = squeeze(nanvar(ts,[],1)); % var_matrix will be epochs x channels
var_grand_median = nanmedian(var_matrix(:)); % grand median

bad_epochs = var_matrix < thresh(1) * var_grand_median | ...
    var_matrix > thresh(2) * var_grand_median;

return




