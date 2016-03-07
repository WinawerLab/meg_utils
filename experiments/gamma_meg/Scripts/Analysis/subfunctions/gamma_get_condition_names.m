function [condition_names, baseline_condition] = gamma_get_condition_names(session_number)
% [condition_names, baseline_condition] = ...
%       gamma_get_condition_names(session_number)
%
% Condition names for MEG Gamma experiment
%
% Inputs: 
%   session_number: integer, for example, the session number for this
%                             folder is 4: 04_Gamma_7_23_2014_subj013
% Outputs:
%  condition_names: cell array of condition names
%
% Example: condition_names = gamma_get_condition_names(5)
%

% TODO fix for session_number <= 3

condition_names               = {   ...
    'White Noise' ...
    'Binary White Noise' ...
    'Pink Noise' ...
    'Brown Noise' ...
    'Gratings(0.36 cpd)' ...
    'Gratings(0.73 cpd)' ...
    'Gratings(1.45 cpd)' ...
    'Gratings(2.90 cpd)' ...
    'Plaid'...
    'Blank'};

if session_number >= 10
    condition_names{3} = 'Binary Pink Noise';
    condition_names{4} = 'Binary Brown Noise';
end

baseline_condition = find(~cellfun(@isempty,strfind(condition_names, 'Blank')));

return
