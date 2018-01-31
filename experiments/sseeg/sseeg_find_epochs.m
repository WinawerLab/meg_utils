function onsets = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run,...
    epochs_per_block)
%% Find the epoch onsets in an SSEEG experiment. 
%
% Find all of the reversals from 0-->1 or 1-->0 in ev_ts. Then, from the
% start index of the run, we make every 12th reversal an epoch timepoint,
% since there are two reversals per DIN event, and six DIN events per epoch.
%
% Inputs: 
%   ev_ts:            a cell (nr_runs x 1) containing arrays of ones and 
%                     zeros (time x DIN value), with length equal to length 
%                     of run (including before and after trial). 
%   images_per_block: number of image reversals occurred in one block 
%   blocks_per_run:   number of blocks occured during one run
%   epochs_per_block: number of epochs in one block 
%
% Outputs:
%   onsets:           a cell (nr_runs x 1) containing arrays of ms timepoints
%                     (framenumber) indicating the onset of every epoch in a 
%                     given run 
%%
nr_runs             = numel(ev_ts);
reversals_per_run   = images_per_block * blocks_per_run;
reversal_inds       = cell(1,nr_runs);
onsets              = cell(1,nr_runs);

for ii = 1:nr_runs
    reversal_inds{ii} = find(diff(ev_ts{ii}))+1;
    onsets{ii}        = reversal_inds{ii}(1:images_per_block/epochs_per_block:reversals_per_run);
end

return