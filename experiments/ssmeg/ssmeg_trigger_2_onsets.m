function onsets = ssmeg_trigger_2_onsets(triggers, which_subject,which_data)
% Find epoch onsets from triggers
%
% For the ssmeg experiment, there were twelve triggers sent per second, and
% we epoch in one second bins
if notDefined('which_data'); which_data = 'meg'; end

switch which_data
    case 'eye'
        inds                = triggers(:,2);
        keep_inds           = inds(1:12:end);
        which_subject       = 7; % If eye data, we need to insert blanks slightly different, 
                                 % and then do nothing like for subject 7
                                 % in meg data.
        
        blanks = [];
        last_onperiod_onset_frames = keep_inds(6:6:end);
        for jj = 1:length(last_onperiod_onset_frames)
            for ii = 1:6
                blanks = [blanks; ii*1000 + last_onperiod_onset_frames(jj)];
            end
        end
        onsets = sort([keep_inds;blanks]);
        
        
    case {'eeg','meg'}
        inds                = find(triggers);
        keep_inds           = inds(1:12:end);
        onsets              = zeros(size(triggers));
        onsets(keep_inds)   = triggers(keep_inds);
end

switch which_subject
    case {7, 8}
        % Done. Do nothing
    case {1,2,3,4,6}
        % Insert blanks every 1000 ms following last on epoch of each block
        last_onperiod_onset_frames = keep_inds(6:6:end);
        for ii = 1:6
            onsets(ii*1000 + last_onperiod_onset_frames) = 3;
        end
        
end

return
