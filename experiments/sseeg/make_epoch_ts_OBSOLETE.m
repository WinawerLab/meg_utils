function epoch_ts = make_epoch_ts(order, nr_runs, ev_ts, epoch_starts)
% Create a timeseries that marks the start of every epoch, as opposed to
% every trigger. Additionally, meg_make_epochs requires a ts with
% stimulus condition numbers (e.g. 1, 3, 5, and 7 for full visual field,
% blank, left field, right field in the SSEEG experiment), instead of just
% using ones to mark events, as we did in ev_ts. 
%   
%   INPUTS: 
% order         : vector representing the stimulus condition sequence.
%                   (e.g. [1 3 1 3 5 3 5 3 7 3 7 3])
% nr_runs       : number of runs in a session 
% ev_ts         : A cell with an event timeseries for each run, padded by
%                   the refresh rate to indicate on and off timeperiods        
% epoch_starts  : A cell with a vector for each run. Each vector contains
%                   the timepoints in seconds of the start of each
%                   epoch, for all epochs, on and off. 

nr_trials = size(order,2);
epochs_per_block = nr_trials/2;

for ii = 1:nr_runs
    epoch_ts{ii} = zeros(1,length(ev_ts{ii})-1);
    for jj = 1:nr_trials
        epoch_ts{ii}(epoch_starts{ii}((jj-1)*epochs_per_block+1:jj*epochs_per_block)) = order(1,jj);
    end
end

return





for ii = 1:nr_runs
    epoch_ts{ii} = zeros(1,length(ev_ts{ii})-1);
    for ll = 1:length(order_long)
        epoch_ts{ii}(epoch_starts{ii}(ll)) = order_long{ii}(ll);
    end
end


for ii = 1:nr_runs
    figure; plot(epoch_ts{ii})
end
