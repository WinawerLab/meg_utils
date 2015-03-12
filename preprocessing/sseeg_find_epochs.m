function [on_epoch_times, off_epoch_times] = sseeg_find_epochs(ev_ts, t, start_inds, ...
    num_time_points, trigs_per_block, blocks_per_run, plot_figures)

% Function to find the timepoints of all 'on' and 'off' epochs in an SSEEG
% experiment. 

%% Find on trigger times
% Find all of the reversals from 0-->1 or 1-->0 in ev_ts. Then, from the
% start index of the run, we make 12th reversal an epoch timepoint, since
% there are two reversals per DIN, and six DIN's per epoch. A couple hacks,
% start_run_ind{ii} is shifted by one from reversal_inds. Additionally, we
% might want a better way of defining the value '431', it is 432-1, 432
% being 72 times 6 (reversals per block times blocks
% per run).
nr_runs         = size(start_inds,2);
trigs_per_run   = trigs_per_block * blocks_per_run;

reversal_inds   = cell(1,nr_runs);
on_epoch_times  = cell(1,nr_runs);
off_starts      = cell(1,nr_runs);

for ii = 1:nr_runs
    reversal_inds{ii} = find(diff(ev_ts{ii}==1));
    on_epoch_times{ii} = t{ii}(reversal_inds{ii}(1:12:(trigs_per_run*2)-1));
    off_starts{ii} = (reversal_inds{ii}(trigs_per_block*2:trigs_per_block*2:end));
end

off_inds_temp = zeros(6,6);
off_epoch_times = cell(1,nr_runs);
for ii = 1:nr_runs
        for jj = 1:blocks_per_run
            off_inds_temp(jj,:) = off_starts{ii}(jj):num_time_points:off_starts{ii}(jj)...
                +((trigs_per_block/blocks_per_run)-1)*num_time_points;
        end
    off_epoch_times{ii} = t{ii}(reshape(off_inds_temp, [1,trigs_per_block]));
end

%%

if plot_figures
    for ii = 1:nr_runs
       figure; plot(t{ii}, ev_ts{ii});
       hold on; plot(on_epoch_times{ii}, 1, 'rx');
       hold on; plot(off_epoch_times{ii}, 1, 'bx');
    end
end

return

