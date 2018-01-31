function [ab, log_frequency, bbPoly, frequencies] = sseeg_broadband(num_epoch_time_pts, amps)
% The ssmeg_broadband function reshapes the frequency amplitude data and fits one line for
% all trials (for the on flicker and off periods). The frequencies with the 
% stimulus locked peaks will be removed before calculating the broadband
% signal between on and off periods.
% 
% INPUTS:
% num_epoch_time_pts           = Number of timepoints in epochs (i.e. 167)
% amps_*_(full, left, right)   = Absolute FFT data for the 6 different
%                                conditions
% OUTPUTS
% ab                           = matrix of broadband amplitudes, num conditions x num channels
% log_frequency                = Array of the used frequencies in log space for calculating
%                               the polynomial fits.
% bbPoly                       = Array containing all the polynomial line
%                                fits per channel (intercept and slope)
% frequencies                  = The used frequencies for calculating the 
%                                broadband, needed for plotting purposes.
%                                
                            
                            
%% Check for code
if isempty(which('eegSummarizeSpectra')),
    addpath(genpath('/Volumes/server/Projects/EEG/Code/')); 
end

if isempty(which('ssmeg_GetSLandABfrequencies')),
    addpath(genpath('/Volumes/server/Projects/MEG/SSMEG/Code/'));
end

num_conditions = size(amps,4);

%% Plot the usual visual channels for full on and off periods

figures = false;
freq = (0:num_epoch_time_pts-1)/(num_epoch_time_pts/1000);

if figures
    
    chan=1;
    figure; clf; hold all;
    for ii = 1:num_conditions
        plot(freq, nanmedian(amps(:,:,chan,ii),2));
        % hold on;
    end
    get_MEG_axes('True');
    xlim([10 150]);
    title(sprintf('Full on and blank periods for channel %d', chan))


    chan=14;
    figure; clf; hold all;
    for ii = 1:num_conditions        
        plot(freq, nanmedian(amps(:,:,chan,ii),2));
        %    hold on;
    end
    get_MEG_axes('True');
    xlim([10 150]);
    title(sprintf('Full on and blank periods for channel %d', chan))
end

%% Define variables
    
conditionNumbers = 1:num_conditions; 


num_channels = 128;
num_freqruencies = numel(freq);


%% Mean over epochs for both conditions

data = zeros(num_freqruencies, num_channels, num_conditions);
for ii = 1:num_conditions
    data(:,:,ii) = squeeze(nanmedian(amps(:,:,1:num_channels,ii),2));
end

%% Reshape to use it in the MEG code

% We want condition x frequency x channels
amps = permute(data, [3 1 2]);


%% Compute broadband

max_time_epoch     = 1000;
freq_for_one_epoch = freq;

T = 1;
f = freq;

fmax = 150; 
% Ph = NaN(size(amps));

% Get the stimulus-locked and broadband frequencies
frequencies = ssmeg_GetSLandABfrequencies(f(f<fmax), T, 12);

%Make empty arrays
ab = zeros(num_conditions, num_channels);  % asynchronous broadband amplitude
bbPoly = zeros(2, num_channels); % broadband polynomial fit (log log space)

% % for debugging: set the broadband frequencies to be the same as the
%       stimulus locked
% frequencies.ab = (1:4) * frequencies.sl ;
% frequencies.ab_i = frequencies.sl_i:frequencies.sl:49;

% loop over channels
% chan = 1;
for ii = 1:128
    
    % fit the spectral data with a line in log-log space
    
    % Take the log of the frequencies and squared amplitude to get log-log space
    log_power     = log(amps(:,frequencies.ab_i, ii).^2);
    log_frequency = repmat(log(frequencies.ab), size(log_power, 1), 1);
    
    % fit one line to all the concatenated trials to get a single slope
    p = polyfit(log_frequency(:), log_power(:), 1);
    
    % evaluate the fitted broadband line at the stimulus-locked frequency
    ssa = polyval(p, log(frequencies.sl));
    
    % calculate the residuals from the linear prediction
    residual = bsxfun(@minus, log_power, polyval(p, log(frequencies.ab)));
    
    % take the mean of the residuals from each condition across frequencies and add the intercept
    % from the linear fit to get the  intercept for that condition
    intercepts = nanmean(residual,2) + ssa;
    
    % broadband based on linear fit (we exponentiate to have units of power,
    % rather than log power)
    ab(:,ii) = exp(intercepts);
    
    % store the polynomial fits
    bbPoly(:,ii) = p;
    
end
% 
% sl      = squeeze(amps_full(:,frequencies.sl_i,:));
% slnoise = squeeze(mean(amps_full(:,[frequencies.sl_i-1 frequencies.sl_i+1],:),2));
% slSNR   = sl ./ slnoise;
% 
% figure;
% plot(log_frequency',polyval(p, log_frequency)'+ab(1,chan), 'b--'); hold on;
% plot(log_frequency',polyval(p, log_frequency)'+ab(2,chan), 'r--'); hold on;
% plot(log_frequency(1,:)', log(nanmedian(amps_on_full(frequencies.ab_i,:,chan),2).^2), 'b'); hold on;
% plot(log_frequency(2,:)', log(nanmedian(amps_off_full(frequencies.ab_i,:,chan),2).^2), 'r'); hold on;
% legend('Flicker fitted line','Blank fitted line', 'Full flicker periods', 'Blank periods' );
% 
% 
% set(gca , 'YScale', 'log', 'XScale', 'log');



% get_MEG_axes('True') % GIVE LOG-LOG axes






    