function onsets = sseeg_find_epochs(ev_ts, images_per_condition, conditions_per_run,...
    epochs_per_condition)
% Find the epoch onsets in an SSEEG experiment. 

% Find all of the reversals from 0-->1 or 1-->0 in ev_ts. Then, from the
% start index of the run, we make ev    ery 12th reversal an epoch timepoint,
% since there are two reversals per DIN, and six DIN's per epoch.
nr_runs             = numel(ev_ts);
reversals_per_run   = images_per_condition * conditions_per_run;
reversal_inds       = cell(1,nr_runs);
onsets              = cell(1,nr_runs);

% off_starts      = cell(1,nr_runs);

for ii = 1:nr_runs
    reversal_inds{ii} = find(diff(ev_ts{ii}))+1;
    onsets{ii}        = reversal_inds{ii}(1:images_per_condition/epochs_per_condition:reversals_per_run);
end

%   this line used to be in the loop but is no longer needed as per below message   
%   off_starts{ii}    = reversal_inds{ii}(trigs_per_block:trigs_per_block:end);

%% probably no longer needed, since we now have a flicker during the off periods
% 
% epoch_length     = median(diff(onsets{1}));
% epochs_per_block = 6;
% 
% for this_run = 1:nr_runs    
%     for this_epoch = 2:epochs_per_block
%         off_starts{this_run}(this_epoch,:) = off_starts{this_run}(this_epoch-1,:)+epoch_length;
%     end
%     off_starts{this_run} = off_starts{this_run}(:)';
%     onsets{this_run} = sort([onsets{this_run}  off_starts{this_run}]);
% end



return

