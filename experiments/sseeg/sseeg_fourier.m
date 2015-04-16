function [amps_on_full,amps_on_right,amps_on_left, amps_off_full,amps_off_right ...
    ,amps_off_left] = sseeg_fourier(t, num_epoch_time_pts, ts, conditions, ...
   num_epochs, data_channels, channels_to_plot, produce_figures)

%% Function description
% This function will calculate the FFT of the different input timeseries.
% And will make some nice plots of the amplitude data.
%
% INPUTS:
%
% t:                    Relative time (0 = start of first trigger)
% num_epoch_time_pts:   Number of timepoints in one epoch (i.e. 167)
% ts:                   Timeseries (1000 x 180 x 192) for flicker
%                       periods (6 different conditions concatenated)
% conditions:           Vector of conditions concatenated across all runs
% num_epochs:           Number of epochs per condition
% channels_to_plot:     Channels to be plotted 

% OUTPUTS:
% amps_*_*:             Absolute fourier transformed amplitudes for every
%                       condition.

% add path for the function 'get_MEG_axes' which formats the plots properly
addpath(genpath('/Volumes/server/Projects/MEG/SSMEG/code/'));


freq                = (0:num_epoch_time_pts-1)/(num_epoch_time_pts/1000); 
off_conditions      = find(conditions ==3);
a                   = size(off_conditions,2)/3;

% find and extract the six different conditions from the timeseries
off_full  = ts(:,off_conditions(1:a), data_channels);
off_right = ts(:,off_conditions(a+1:2*a), data_channels);
off_left  = ts(:,off_conditions((2*a)+1:3*a), data_channels);
on_full   = ts(:, find(conditions == 1), :);   
on_right  = ts(:, find(conditions == 5), :);
on_left   = ts(:, find(conditions == 7), :);   

% calculate fast fourier transform
ft_on_epoched_full   = fft(on_full) / length(t)*2;   
ft_on_epoched_right  = fft(on_right) / length(t)*2;
ft_on_epoched_left   = fft(on_left) / length(t)*2; 
ft_off_epoched_full  = fft(off_full)/ length(t)*2; 
ft_off_epoched_right = fft(off_right)/ length(t)*2; 
ft_off_epoched_left  = fft(off_left)/ length(t)*2; 

% take absolute values
amps_on_full   =  abs(ft_on_epoched_full);    clear ft_on_epoched_full;
amps_on_right  =  abs(ft_on_epoched_right);   clear ft_on_epoched_right;
amps_on_left   =  abs(ft_on_epoched_left);    clear ft_on_epoched_left;
amps_off_full  =  abs(ft_off_epoched_full);   clear ft_off_epoched_full;
amps_off_right =  abs(ft_off_epoched_right);  clear ft_off_epoched_right;
amps_off_left  =  abs(ft_off_epoched_left);   clear ft_off_epoched_left;

if produce_figures;
    
% Blank periods
    figure(203);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    plot(freq, squeeze(nanmedian(amps_off_full(:,:,channels_to_plot), 2)), 'o-')
    xlim([3 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Blank periods')
    get_MEG_axes('True');

    %% All three blank periods (left, right, full)

    figure(204);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20)
    hold on;
    a = squeeze(nanmedian(amps_off_full(:,:,channels_to_plot), 2));
    b = squeeze(nanmedian(amps_off_left(:,:,channels_to_plot), 2));
    c = squeeze(nanmedian(amps_off_right(:,:,channels_to_plot), 2));
    plot(freq, a, freq, b, freq, c, 'LIneWidth', 2)
    xlim([3 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Blank periods')
    get_MEG_axes('True');

    %% Full, left, right field on

    figure(206);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    plot(freq, squeeze(nanmedian(amps_on_full(:,:,channels_to_plot), 2)), 'LineWidth', 2)
    xlim([3 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker periods for full field stimulus')
    get_MEG_axes('True');

    % Right half on
    figure(207);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    plot(freq, squeeze(nanmedian(amps_on_right(:,:,channels_to_plot), 2)), 'LineWidth', 2)
    xlim([3 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker periods for presenting stimulus right')
    get_MEG_axes('True');

    % Left half on
    figure(208);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    plot(freq, squeeze(nanmedian(amps_on_left(:,:,channels_to_plot), 2)), 'LineWidth', 2)
    xlim([3 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker periods for presenting stimulus left')
    get_MEG_axes('True');


    %% Difference plots

    % Difference Full field and blank
    figure(209);
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    hold on;
    plot(freq, squeeze(nanmedian(amps_on_full(:,:,channels_to_plot), 2)) ./ squeeze(nanmedian(amps_off_full(:,:,channels_to_plot), 2)))
    xlim([3 110])
    ylim([-.5 4.5])
    yl = get(gca, 'YLim');
    for ii = 1:10; plot(ii*freq(2)*[12 12], yl, 'k-'); end

    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker full field minus blank')

    % Difference right field and blank
    figure(210);
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    hold on;
    plot(freq, squeeze(nanmedian(amps_on_right(:,:,channels_to_plot), 2)) ./ squeeze(nanmedian(amps_off_right(:,:,channels_to_plot), 2)))
    xlim([3 110])
    ylim([-.5 4.5])
    yl = get(gca, 'YLim');
    for ii = 1:10; plot(ii*freq(2)*[12 12], yl, 'k-'); end

    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker right field minus blank')

    % Difference left field and blank
    figure(211);
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    hold on;
    plot(freq, squeeze(nanmedian(amps_on_left(:,:,channels_to_plot), 2)) ./ squeeze(nanmedian(amps_off_left(:,:,channels_to_plot), 2)))
    xlim([3 110])
    ylim([-.5 4.5])
    yl = get(gca, 'YLim');
    for ii = 1:10; plot(ii*freq(2)*[12 12], yl, 'k-'); end

    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    title('Flicker left field minus blank')

    end
    
end