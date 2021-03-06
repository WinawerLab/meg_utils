function [start_ind] = eeg_find_trial_begin(init, event_ts)
%

% Cross correlate the flicker sequence and the event time series to find
% events which correspond to the photodiode pre-scan sequence
[c, lags] = xcov(double(event_ts), double(init));

% The peak of the cross-correlation indicates the alignment
[~, tmp] = max(c);
start_flicker_ind = lags(tmp);
end_flicker_ind   = start_flicker_ind + length(init);

% Zero out all events up to the end of the photodiode pre-scan sequence
event_ts(1:end_flicker_ind+20) = 0;
start_ind = find(event_ts,1);

% return

% debug
if nargout == 0,
    figure
    plot(event_ts, 'r'); hold on
        axis([min(axis) max(axis) 0 1.5]);
    inds_to_plot = start_flicker_ind+(0:length(init)-1);
    plot(inds_to_plot, event_ts(inds_to_plot), 'g')
    plot(inds_to_plot, init, 'k--')
    
%     xlim([inds_to_plot(1)-1000 inds_to_plot(end)+1000])
    plot(start_ind * [1 1], [0 1], 'b--', 'LineWidth',3)
end