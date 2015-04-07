%% EEG script to plot values on a topoplot

data_channels = 1:128;
ts  = ts_one(:,:,data_channels);
fs  = 1000; % sample rate
T   = .995; %s
slF = 12; %Hz

% Compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);


% Compute fft and mean across epochs per channel
% spectral_data = abs(fft(ts))/length(t)*2;
% spectral_data_mean = squeeze(nanmean(spectral_data,2));
% spectral_data_std  = squeeze(nanstd(spectral_data_boots, [], 2));
% spectral_data_snr  = spectral_data_mean./spectral_data_std;

%[time x epochs x epochs] --> [channels x time x epochs]
ts = permute(ts, [3 1 2]);
% Calculate fft and SL power the GLM Denoise way
freq = megGetSLandABfrequencies(f, T, slF);
sl   = getstimlocked(ts,freq);

mean_sl = mean(sl);
clims   = [0 25];
ttl     = 'Stimulus locked, all trials';

eegPlotMap2(mean_sl,clims,[],[],ttl,[])


