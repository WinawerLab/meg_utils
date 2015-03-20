function epoch_ts = make_epoch_ts(order, nr_runs, ev_ts, epoch_starts)
% Create a timeseries that marks the start of every epoch, as opposed to
% every trigger. Additionally, meg_make_epochs requires a ts with
% stimulus condition numbers (e.g. 1, 3, 5, and 7 for full visual field,
% blank, left field, right field in the SSEEG experiment), instead of just
% using ones to mark events, as we did in ev_ts. 
%   
%   INPUTS: 
% order             : a vector representing the stimulus condition sequence.
%       (e.g. [1 3 1 3 5 3 5 3 7 3 7 3])

nr_trials = size(order,2);
epochs_per_block = nr_trials/2;

for ii = 1:nr_runs
    epoch_ts{ii} = zeros(1,length(ev_ts{ii})-1);
    for jj = 1:nr_trials
        epoch_ts{ii}(epoch_starts{ii}((jj-1)*epochs_per_block+1:jj*epochs_per_block)) = order(1,jj);
    end
end

return
