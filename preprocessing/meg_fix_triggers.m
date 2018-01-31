function trigger = meg_fix_triggers(ts)
% Convert an analog signal representing triggers in binary, and return a
% vector of integer trigger values over time.
%
%   trigger = meg_fix_triggers(ts)
%
% Inputs:
%   ts: matrix of trigger values (time points by channel). Each column is
%       an analog representastion of a binary value (0 or 1). We thus convert
%       the columns to binary and then combine the values in a row from
%       binary to a base 10 integer. 
%
% Output: vector of base 10 trigger values over time (time points x 1)
%
% Example:
%       ts = [0 0 0; 0 0 0; 1 0 0; 0 0 0; 1 1 1];
%       trigger = meg_fix_triggers(ts)
%       This should return trigger = [0 0 1 0 7];

% error('Not yet implemented')

% exclude outliers
md = median(ts(:));
ts(ts > 5*md) = NaN;

 
% rescale to [0 1]
ts = ts - min(ts(:));
ts = ts / max(ts(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(ts(:))) == 0, trigger_is_high = true; 
else                     trigger_is_high = false; end 

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, ts = 1 - ts; end

% threshold to binarize from analog recording 
ts = ts > 0.1;
% differentiate to isolate trigger onsets, and pad to preserve length
ts = padarray(diff(ts), [1 0], 0, 'post');

% rectify 
ts(ts<0) = 0;

% mark the samples as 1 when any trigger channel signaled 
any_trigger      = sum(ts,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);

% check whether 2 triggers signalled in an impossibly short time period
% (indicating that they were probably simultanuous, but slighly misaligned)
triggers_that_are_too_close = diff(any_trigger_inds) < 10;

time_points_that_are_too_close = any_trigger_inds(triggers_that_are_too_close);

% Issue: most of the above indexes return all zeros, check on that 
% Solution: because of padding, the triggers are present -1 or -2 time
% points before the the time points listed in bad_time_points

%% Realign bad triggers
%   We find any case where two trigger channel onsets were within 2
%   samples, and we re-align them to the first sample.
for ii=1:length(time_points_that_are_too_close)
    ts(time_points_that_are_too_close(ii),:) = sum(ts(time_points_that_are_too_close(ii)+[-2:2],:));
    ts(time_points_that_are_too_close(ii) + [-2 -1 1 2],:) = 0;
end

% now check that there are no time points that are too close
any_trigger                 = sum(ts,2) > 0;
any_trigger_inds            = find(any_trigger);
triggers_that_are_too_close = diff(any_trigger_inds) < 10;
assert(sum(triggers_that_are_too_close) == 0)

%% convert binary matrix into base 10 vector
trigger = ts * 2.^(0:size(ts,2)-1)';


