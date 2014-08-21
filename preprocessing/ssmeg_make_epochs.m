function [ts_on,ts_off,num_epochs] = ssmeg_make_epochs(ts, num_epoch_time_pts,epochs_on,epochs_off)

% This function makes epochs (timepoints x nr of epochs x nr of channels),
% by the use of the indices made by ssmeg_fix_triggers.

% INPUTS:
% ts:                   timeseries
% num_epoch_time_pts:   Number of timepoints per single epoch (i.e. 167)
% epochs_on:            Array with start indices of the three 
%                       flickering conditions (Full field, left and right)
% epochs_off:           Array with start indices of the three 
%                       blank periods (after a full field, left and right flicker period)

% OUTPUTS:
% ts_on:                Timeseries for 

%% Define indices per condition

% ON periods
epoch_start_full_on_ind  = epochs_on(:,1);
epoch_start_right_on_ind = epochs_on(:,2);
epoch_start_left_on_ind  = epochs_on(:,3);

% OFF periods
epoch_start_full_off_ind  = epochs_off(:,1);
epoch_start_right_off_ind = epochs_off(:,2);
epoch_start_left_off_ind  = epochs_off(:,3);

%% Prepare for making epochs
num_channels = size(ts, 2);
num_epochs   = length(epoch_start_full_on_ind); % The same for every condition (i.e. 72 time pts)

% Make new arrays
ts_on_full_epoched   = zeros(num_epoch_time_pts, num_epochs, num_channels); % Time x Epochs x Channels
ts_on_right_epoched  = zeros(num_epoch_time_pts, num_epochs, num_channels); 
ts_on_left_epoched   = zeros(num_epoch_time_pts, num_epochs, num_channels);
ts_off_full_epoched  = zeros(num_epoch_time_pts, num_epochs, num_channels); 
ts_off_right_epoched = zeros(num_epoch_time_pts, num_epochs, num_channels); 
ts_off_left_epoched  = zeros(num_epoch_time_pts, num_epochs, num_channels);


%% Loop over channels and epochs of every condition

for channel = 1:num_channels
    
   
    for epoch = 1:num_epochs
        
        % **ON PERIODS** 
        
        %Full field
        ts_on_full_epoched(:, epoch, channel) = ...
            ts(epoch_start_full_on_ind(epoch) + (0:num_epoch_time_pts-1), channel);
        
        % Right field
        ts_on_right_epoched(:, epoch, channel) = ...
            ts(epoch_start_right_on_ind(epoch) + (0:num_epoch_time_pts-1), channel);
        
        % Left field
        ts_on_left_epoched(:, epoch, channel) = ...
            ts(epoch_start_left_on_ind(epoch) + (0:num_epoch_time_pts-1), channel);       
        
        % **OFF PERIODS** 
        
        %Full field
        ts_off_full_epoched(:, epoch, channel) = ...
            ts(epoch_start_full_off_ind(epoch) + (0:num_epoch_time_pts-1), channel);
        
        % Right field
        ts_off_right_epoched(:, epoch, channel) = ...
            ts(epoch_start_right_off_ind(epoch) + (0:num_epoch_time_pts-1), channel);
        
        % Left field
        ts_off_left_epoched(:, epoch, channel) = ...
            ts(epoch_start_left_off_ind(epoch) + (0:num_epoch_time_pts-1), channel);                             
    end

end

ts_on  = [ts_on_full_epoched,ts_on_right_epoched,ts_on_left_epoched];
ts_off = [ts_off_full_epoched,ts_off_right_epoched,ts_off_left_epoched];


clear ts