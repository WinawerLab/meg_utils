function [bad_epochs_on,bad_epochs_off,ts_bad_clean_on,ts_bad_clean_off] = ...
    ssmeg_bad_epochs(t,num_epoch_time_pts, ts_on,ts_off, num_epochs, produce_figures)

% This function will take every epoch and compares the variance of the
% signal within this epoch with a certain threshold (say, 10 times the
% median of the signal of one channel). All epochs with a variance larger
% than the threshold, will be defined as a bad epoch (1) in a matrix. One 
% can use this matrix to plot the difference between with/without bad
% epochs.

% INPUTS
% t:                    Relative time
% num_epoch_time_pts:   Number of timepoints in epoch (i.e. 167)
% ts_for_epochs_on:     Timeseries for flicker on periods (3 conditions
%                       concatenated)
% ts_for_epochs_off:    Same, but then for off (blank) periods (3 conditions
%                       concatenated).

% OUTPUTS 
% bad_epochs_o*:        Matrix with all the bad_epochs in on/off periods
% ts_bad_clean_*:       New timeseries with NaN's for the bad epochs

%%

%% Define timeseries again
ts_clean_on_full   = ts_on(:,1:num_epochs,:);
ts_clean_on_right  = ts_on(:,(num_epochs+1):2*num_epochs,:);
ts_clean_on_left   = ts_on(:,(2*num_epochs)+1:3*num_epochs,:);

ts_clean_off_full   = ts_off(:,1:num_epochs,:);
ts_clean_off_right  = ts_off(:,(num_epochs+1):2*num_epochs,:);
ts_clean_off_left   = ts_off(:,(2*num_epochs)+1:3*num_epochs,:);



%% Plot two visual channels
if produce_figures

    chan_1 = 1;
    figure(106); plot(squeeze(mean(ts_clean_on_full(:,:,chan_1),2)));
    xlabel('Time (ms)')
    ylabel('Amplitude (Tesla)')
    title(sprintf('Timeseries of channel nr %d', chan_1))

    chan_14 = 14;
    figure(107); plot(squeeze(mean(ts_clean_off_full(:,:,chan_14),2)));
    xlabel('Time (ms)')
    ylabel('Amplitude (Tesla)')
    title(sprintf('Timeseries of channel nr %d', chan_14))

end

%% Make matrix to remove bad epochs for ON-data

% epochs X channel
bad_epochs_on_full_clean=zeros(size(ts_clean_on_full,2),size(ts_clean_on_full,3));
bad_epochs_off_full_clean=zeros(size(ts_clean_off_full,2),size(ts_clean_off_full,3));

bad_epochs_on_left_clean=zeros(size(ts_clean_on_left,2),size(ts_clean_on_left,3));
bad_epochs_off_left_clean=zeros(size(ts_clean_off_left,2),size(ts_clean_off_left,3));

bad_epochs_on_right_clean=zeros(size(ts_clean_on_right,2),size(ts_clean_on_right,3));
bad_epochs_off_right_clean=zeros(size(ts_clean_off_right,2),size(ts_clean_off_right,3));

remove_fun = @(x) x>(10*median(x)) | x < (.1 * median(x));

% check epochs
for chan=1:size(ts_clean_on_full,3)%channels   
    % FULL ON
    a=nanvar(squeeze(ts_clean_on_full(:,:,chan)),[],1);
    bad_epochs_on_full_clean(remove_fun(a),chan)=1;
        
    % FULL OFF
    a=nanvar(squeeze(ts_clean_off_full(:,:,chan)),[],1);
    bad_epochs_off_full_clean(remove_fun(a),chan)=1;
        
    % LEFT ON
    a=nanvar(squeeze(ts_clean_on_left(:,:,chan)),[],1);
    bad_epochs_on_left_clean(remove_fun(a),chan)=1;
    
    % LEFT OFF
    a=nanvar(squeeze(ts_clean_off_left(:,:,chan)),[],1);
    bad_epochs_off_left_clean(remove_fun(a),chan)=1;
    
    % RIGHT ON
    a=nanvar(squeeze(ts_clean_on_right(:,:,chan)),[],1);
    bad_epochs_on_right_clean(remove_fun(a),chan)=1;
    
    % RIGHT OFF
    a=nanvar(squeeze(ts_clean_off_right(:,:,chan)),[],1);
    bad_epochs_off_right_clean(remove_fun(a),chan)=1;
    
end

% check channels (I would do this in a different way)
% for k=1:size(ts_on_epoched,1)%epochs
%     a=var(squeeze(ts_on_epoched(k,:,:)),[],1);
%     bad_epochs_on(k,a>(10*median(a)))=1;
% end

%% plot bad epoch matrix for on and off data
if produce_figures

    figure
    imagesc(bad_epochs_on_full_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');

    figure
    imagesc(bad_epochs_off_full_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');
    
    figure
    imagesc(bad_epochs_on_left_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');
    
    figure
    imagesc(bad_epochs_off_left_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');
    
    figure
    imagesc(bad_epochs_on_right_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');
    
   
    figure
    imagesc(bad_epochs_off_right_clean')
    xlabel('Epochs');
    ylabel('Channels');
    title('Removed bad epochs');
    legend('Removed = Red / Not removed = Blue');
    
    
end


%% Plot two visual channels again

if produce_figures

    figure(107)
    chan_1 = 98;

    temp_data=ts_clean_on_full(:,:,chan_1);
    plot(nanmean(temp_data(:,bad_epochs_on_full_clean(:,chan_1)==0),2), 'r'); hold on;
    plot(nanmean(ts_clean_on_full(:,:,chan_1),2),'b');
    xlabel('Time (ms)');
    ylabel('Amplitude (Tesla)');
    title(sprintf('Timeseries with removed bad epochs of channel nr %d', chan_1));
    legend('Removed','Not removed');

    figure(109)
    chan_14 = 14;

    temp_data=ts_clean_on_full(:,:,chan_14);
    plot(nanmean(temp_data(:,bad_epochs_on_full_clean(:,14)==0),2), 'r'); hold on;
    plot(nanmean(ts_clean_on_full(:,:,14),2),'b');
    xlabel('Time (ms)');
    ylabel('Amplitude (Tesla)');
    title(sprintf('Timeseries with removed bad epochs of channel nr %d', chan_14));
    legend('Removed','Not removed');

end
%% Define clean data with removed epochs as a new variable, with for removed epochs NaN's

 
% Loop over channels & epochs
for chan = 1:157

    % FULL
    for epoch = 1:length(bad_epochs_on_full_clean(:,chan))
        % If a bad epoch array is defined as '1', it's a bad epoch
        if bad_epochs_on_full_clean(epoch,chan) == 1
            % Turn this epoch into NaN's
            ts_clean_on_full(:,epoch,chan) = NaN;
        end

        % Same for off periods
        if bad_epochs_off_full_clean(epoch,chan) == 1
            ts_clean_off_full(:,epoch,chan) = NaN;
        end
    end

    % LEFT
    for epoch = 1:length(bad_epochs_on_left_clean(:,chan))
        % If a bad epoch array is defined as '1', it's a bad epoch
        if bad_epochs_on_left_clean(epoch,chan) == 1
            % Turn this epoch into NaN's
            ts_clean_on_left(:,epoch,chan) = NaN;
        end

        % Same for off periods
        if bad_epochs_off_left_clean(epoch,chan) == 1
            ts_clean_off_left(:,epoch,chan) = NaN;
        end
    end  


    % RIGHT
    for epoch = 1:length(bad_epochs_on_right_clean(:,chan))
        % If a bad epoch array is defined as '1', it's a bad epoch
        if bad_epochs_on_right_clean(epoch,chan) == 1
            % Turn this epoch into NaN's
            ts_clean_on_right(:,epoch,chan) = NaN;
        end

        % Same for off periods
        if bad_epochs_off_right_clean(epoch,chan) == 1
            ts_clean_off_right(:,epoch,chan) = NaN;
        end
    end  

end

%% Redefine variables
ts_bad_clean_on     = [ts_clean_on_full, ts_clean_on_right, ts_clean_on_left,];
ts_bad_clean_off    = [ts_clean_off_full, ts_clean_off_right, ts_clean_off_left];


bad_epochs_on       = [bad_epochs_on_full_clean,bad_epochs_on_right_clean,bad_epochs_on_left_clean];
bad_epochs_off      = [bad_epochs_off_full_clean,bad_epochs_off_right_clean, bad_epochs_off_left_clean];