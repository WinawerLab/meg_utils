function F = gamma_broadband_fit_loglog(x,P,f,p_exp)

% function for fitting log power log frequnecy spectrum as line + gaussian
%
%     Copyright (C) 2014  D Hermes
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

% P:        data (log10(power))
% f:        log10(frequencies)
% p_exp:    power law exponent of baseline condition
% x(1):     power law exponent after subtracting baseline
% x(2):     gaussian height (log10 power)
% x(3):     gaussian center frequency (log10 Hz)
% sigma:    gaussian standard deviation (log10 Hz)
%
% Residual (F) = ...
%   Data - Prediction
%
% Data: P
%
% Prediction: broadband offset + broadband baseline + gaussian
%   broadband offset: x(1)
%   broadband baseline: -p_exp*f
%   gaussian: normpdf(f,x(3),sigma))

sigma = 0.04;
F= P - (x(1)-p_exp*f + x(2)*.04*sqrt(2*pi)*normpdf(f,x(3),sigma));
% note: .04*sqrt(2*pi) gives an amplitude of 1 to the Gaussian for x(2)=1;
