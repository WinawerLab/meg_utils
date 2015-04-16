function onsets = sseeg_find_epochs(ev_ts, images_per_block, blocks_per_run,...
    epochs_per_block)
% Find the epoch onsets in an SSEEG experiment. 
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

%% INSTRUCTIONS 
%   For a regular run with DIN events during both on and off periods:
%   In Section I, comment out the line starting with "off_starts{ii}" in the first loop, 
%   remove the "/2" from the second line starting with "onsets{ii}" in the 
%   first loop, and comment out everything in Section II. 
%   
%   For a run with DIN events occurring only during on periods:
%   Do the reverse of the above instructions. 

%% Section I

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


%              
%% Section II - Comment this section if DINs are sent during all blocks (on and off)

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

