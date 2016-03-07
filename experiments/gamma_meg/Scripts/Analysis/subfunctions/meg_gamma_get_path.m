function pth = meg_gamma_get_path(session_number)
% pth = meg_gamma_get_path(session_number)
project_pth = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth    = sprintf('%02d_Gamma_*subj*', session_number);

d = dir(fullfile(project_pth, data_pth));

if isempty(d), error('Data path not found'); end
if length(d)>1, error('More than one path found'); end

pth = fullfile(project_pth, d.name);

return