function [baseline_fit,w_broadband,w_gauss,gauss_f,fit_f2] = ...
    gamma_fit_data_localregression_multi(f,f_use4fit,data_base,data_fit)

% function fits broadband + gaussian
% [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
%     gamma_fit_data_localregression_multi(f,f_use4fit,data_base,data_fit);
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
x_base = data_base(f_sel)'; % 1 x num frequencies

% Data for fitting condition of interest
x_in   = data_fit(f_sel,:)'; % num_conditions x num frequencies

% fit baseline with local linear regression, 1 x num_frequencies
baseline_fit = localregression(f_use4fit,log(x_base)',f_use4fit,[],[], []);

x_to_fit = bsxfun(@minus, log(x_in), baseline_fit); % 1 x num_frequencies

my_options=optimset('Display','off','Algorithm','trust-region-reflective');

% F = gamma_broadband_fit_loglog(x,P,f,p_exp)
% Previously (e.g., Hermes 2014), we assumed a Gaussian with sd 0.04 log10
% units. We are  now using natural log. To keep the same bandwidth, we
% multiply sigma by log(10)
sigma = 0.04*log(10); 

num_conditions = size(x_in,1);
idx_gw = 1+(1:num_conditions);
idx_bw = 1+num_conditions+(1:num_conditions);
dummy  = ones(1,num_conditions*2);

[x]=lsqnonlin(@(x) gamma_broadband_fit(x, x_to_fit(:),log(f_use4fit),sigma, idx_gw, idx_bw),...
    [log(50) dummy*0],[log(50) -Inf*dummy ],[log(60) Inf *dummy ],...
    my_options);

gauss_f      = x(1); % gaussian center frequency (log Hz)
w_broadband  = x(idx_bw); % DC offset from baseline (broadband, log power)
w_gauss      = x(idx_gw); % gaussian height (log power)


% fit to data in log-space
fit_f2 = NaN(length(f), num_conditions);

G = (sigma*sqrt(2*pi)*normpdf(log(f_use4fit),gauss_f,sigma))' * w_gauss ;
BBplusG = bsxfun(@plus, G, w_broadband);
BLplusBBplusG = bsxfun(@plus, BBplusG, baseline_fit');
fit_f2(f_sel, :) = BLplusBBplusG;
fit_f2 = fit_f2';
baseline_fit = squeeze(baseline_fit);

function F = gamma_broadband_fit(x,D,f,sigma, idx_gw, idx_bw)
% D:         data (log power)
% f:         frequencies (log Hz)
% x(1):      Gaussian center frequency (log Hz)
% x(idx_bw): DC offset from baseline (broadband, log power, 1xnum_conditions)
% x(idx_gw): gaussian height (log power, 1xnum_conditions)
% sigma:     gaussian standard deviation (log Hz)
% idx_gw:    indices to gaussian height parameters in x
% idx_bw:    indices to broadband parameters in x
%
% Residual (F) = ...
%   Data - Prediction
%
% Data: D
%
% Prediction: offset + gaussian
%   offset: x(1)
%   gaussian: normpdf(f,x(2),sigma))

G = sigma*sqrt(2*pi)*x(idx_gw)'* normpdf(f,x(1),sigma);
P = bsxfun(@plus, G, x(idx_bw)');
F = D - P(:);
% note: sigma*sqrt(2*pi) gives an amplitude of 1 to the Gaussian for x(2)=1;

