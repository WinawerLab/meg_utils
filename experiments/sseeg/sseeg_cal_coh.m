function [coh_full, coh_right, coh_left] = sseeg_cal_coh...
    (amps_on_full,amps_on_right,amps_on_left, plotType, keep_epochs, produce_figures)
% This function will calculate the coherence and plot it in 2D and 3D
% topo plots. The coherence is calculated by taking the signal to noise
% ratio for the stimulus locked peak (signal) and peak of surrounding
% frequencies (noise)
%
% INPUTS:
% amps_*_*:     amplitudes of interest
% meg_files:    filename of the MEG data, you need this for for the header
% plotType: any of '2d', '3d', 'both'

% check inputs
if ~exist('plotType', 'var') || isempty(plotType), plotType = 'both'; end
if ~exist('produce_figures', 'var') || isempty(produce_figures), produce_figures = false; end

%% Define signal and noise frequencies for coherence calculations
%   Noise frequencies are +/- 1 each signal frequency
signal_frequencies = [12 24 36 48];
% signal_frequencies = [30];
noise_frequencies  = sort([signal_frequencies signal_frequencies-1 signal_frequencies+1]);
ss_inds            = signal_frequencies + 1;
noise_inds         = noise_frequencies + 1;


if size(amps_on_full,3) > 128;
    num_chans = 128;
else
    num_chans = size(amps_on_full,3);
end

channelsToPlot = [1:128];

% Full field coherence
signal   = squeeze(sum(amps_on_full(ss_inds,:,channelsToPlot),1));
noise    = squeeze(sum(amps_on_full(noise_inds,:,channelsToPlot),1));
coh_full = nanmean(signal./noise);

% Right field coherence
signal      = squeeze(sum(amps_on_right(ss_inds,:,channelsToPlot),1));
noise       = squeeze(sum(amps_on_right(noise_inds,:,channelsToPlot),1));
coh_right   = nanmean(signal./noise);

% Left field coherence
signal      = squeeze(sum(amps_on_left(ss_inds,:,channelsToPlot),1));
noise       = squeeze(sum(amps_on_left(noise_inds,:,channelsToPlot),1));
coh_left    = nanmean(signal./ noise);

%% Plot coherehnce of stimulus locked signals on mesh

figure(100);
plotOnEgi(coh_full);
title('coherence for full field flicker');

figure(101);
plotOnEgi(coh_right);
title('coherence for right visual flicker');

figure(102);
plotOnEgi(coh_left);
title('coherence for left visual flicker');
