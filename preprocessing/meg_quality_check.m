% MEG Quality check

% What do we want to know?
% 
% Some basics:
%   Recording frequency
%   Recording duration
%   Filters (e.g., high-pass, line noise)
%   Triggers (are they present, and how many? and what values?) (channels 161:168)
%   Was there eye tracking?
%   Was there a photodiode? (channel 191)
%   Units (e.g., Tesla)?
%
% Power spectra of 3 reference channels (158:160)
%   Best fit 1/f?
%   Line noise
%   Other large peaks?
%   Jumps???
%
% How much line noise?
%   the power (relative to what??)
%   the distribution across channels and time
%   which harmonics?
%
% Jumps, artifacts, etc 
%
% Other weird peaks in data
%
% Saturated channels
%   which, when, how many?
%
% Save out images of power spectrum and spectrogram of every channel
%% Look at empty room data
pth = '/Volumes/server/Projects/MEG/Gamma_BR/emptyRoomData/Empty_Room_Gamma Settings_6.2.2016_1135am.sqd';
hdr = ft_read_header(pth);
data = ft_read_data(pth);

