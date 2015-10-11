function [baseline_fit,w_broadband,w_gauss,gauss_f,fit_f2] = ...
    gamma_fit_data_localregression(f,f_use4fit,data_base,data_fit)

% function fits broadband + gaussian
% [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
%     gamma_fit_data_localregression(f,f_use4fit,data_base,data_fit);
%
% input:
%   f:          frequencies in spectral data (Hz)
%   f_use4fit:  frequencies to use for fitting (subset of f, Hz)
%   data_base:  data used to fit spectral shape across conditions without gamma - enter power (not log power)
%   data_fit:   used to fit broadband and gaussian - enter power (not log power)
%
% output (exp weight_pwr weight_gauss gamma_freq fit_f2)
%
%
% Example ??
%


% Indices of the frequencies to use for fitting
f_sel  = ismember(f,f_use4fit);

% Data for fitting baseline (presumably the mean across all conditions)
x_base = data_base(f_sel); % 1 x num frequencies

% Data for fitting condition of interest
x_in   = data_fit(f_sel); % 1 x num frequencies

% fit baseline with local linear regression, 1 x num_frequencies
baseline_fit = localregression(f_use4fit,log(x_base),f_use4fit,[],[], []);
% % no smoothing for baseline
% baseline_fit = log(x_base);

x_to_fit = log(x_in) - baseline_fit; % 1 x num_frequencies

my_options=optimset('Display','off','Algorithm','trust-region-reflective');

% F = gamma_broadband_fit_loglog(x,P,f,p_exp)
% Previously (e.g., Hermes 2014), we assumed a Gaussian with sd 0.04 log10
% units. We are  now using natural log. To keep the same bandwidth, we
% multiply sigma by log(10)
sigma = 0.04*log(10); 

[x]=lsqnonlin(@(x) gamma_broadband_fit(x, x_to_fit,log(f_use4fit),sigma),...
    [0 0 log(50)],[-Inf -Inf log(50)],[Inf Inf log(60)],...
    my_options);

w_broadband  = x(1); % DC offset from baseline (broadband, log power)
w_gauss      = x(2); % gaussian height (log power)
gauss_f      = x(3); % gaussian center frequency (log Hz)

% fit to data in log-space
fit_f2 = NaN(size(f));
fit_f2(f_sel) = w_broadband + baseline_fit + ...
    w_gauss*sigma*sqrt(2*pi)*normpdf(log(f_use4fit),gauss_f,sigma);


function F = gamma_broadband_fit(x,D,f,sigma)
% D:        data log(power)
% f:        log(frequencies)
% x(1):     DC offset from baseline (broadband, log power)
% x(2):     gaussian height (log power)
% x(3):     gaussian center frequency (log Hz)
% sigma:    gaussian standard deviation (log Hz)
%
% Residual (F) = ...
%   Data - Prediction
%
% Data: D
%
% Prediction: offset + gaussian
%   offset: x(1)
%   gaussian: normpdf(f,x(2),sigma))

F = D - (x(1) + x(2)*sigma*sqrt(2*pi)*normpdf(f,x(3),sigma));
% note: sigma*sqrt(2*pi) gives an amplitude of 1 to the Gaussian for x(2)=1;

