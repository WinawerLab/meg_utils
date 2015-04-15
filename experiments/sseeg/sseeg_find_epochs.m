function onsets = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run,...
    epochs_per_block)
% Find the epoch onsets in an SSEEG experiment. 
%
% Find all of the reversals from 0-->1 or 1-->0 in ev_ts. Then, from the
% start index of the run, we make every 12th reversal an epoch timepoint,
% since there are two reversals per DIN event, and six DIN events per epoch.
%
%   INPUTS: 
% ev_ts:            a cell (nr_runs x 1) containing arrays of ones and 
%                   zeros (time x DIN value), with length equal to length 
%                   of run (including before and after trial). 
% images_per_block: number of image reversals occurred in one block 
% blocks_per_run:   number of blocks occured during one run
% epochs_per_block: number of epochs in one block 
%
%   OUTPUTS:
% onsets:           a cell (nr_runs x 1) containing arrays of ms timepoints
%                   (framenumber) indicating the onset of every epoch in a 
%                   given run 

nr_runs             = numel(ev_ts);
reversals_per_run   = images_per_block * blocks_per_run;
reversal_inds       = cell(1,nr_runs);
onsets              = cell(1,nr_runs);

off_starts      = cell(1,nr_runs);

for ii = 1:nr_runs
    reversal_inds{ii} = find(diff(ev_ts{ii}))+1;
    onsets{ii}        = reversal_inds{ii}(1:images_per_block/epochs_per_block:reversals_per_run/2);
    off_starts{ii}    = reversal_inds{ii}(images_per_block:images_per_block:end);
end

%   for a regular run with triggers sent during both on and off periods,
%   take the line that starts with off_starts out of the above loop, take
%   out the '/2' at the end of onsets in above loop to change back to
%   normal, and comment out everything below this message.
%              
%% probably no longer needed, since we now have a flicker during the off periods
% 
epoch_length     = median(diff(onsets{1}));
epochs_per_block = 6;

for this_run = 1:nr_runs    
    for this_epoch = 2:epochs_per_block
        off_starts{this_run}(this_epoch,:) = off_starts{this_run}(this_epoch-1,:)+epoch_length;
    end
    off_starts{this_run} = off_starts{this_run}(:)';
    onsets{this_run} = sort([onsets{this_run}  off_starts{this_run}]);
end



return

