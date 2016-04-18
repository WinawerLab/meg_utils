function parameters = gamma_get_parameters(session_num)
%% parameters = gamma_get_experimental_parameters(session_num)
%  outputs a structure with the experimental parameters corresponding to
%  the given experiment

%% 
parameters = struct('group',[],'condition_names',[],...
    'epoch_range',[],'baseline_condition',[],'ITI',[]);


if session_number < 5
    condition_group = 1;
elseif ismember(session_number, 5:9)
    condition_group = 2;
elseif ismember(session_number, 10:16)
    condition_group = 3;
elseif session_number > 16 && session_number <99
    condition_group = 4;
end


switch condition_group
    case 1
        parameters.group                 = '500ms non-binarized noise';
        parameters.epoch_range           = [0.050 0.549];
        parameters.condition_names       = {   ...
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
        parameters.ITI                   = [];
        
    case 2
        parameters.group                 = '1000ms non-binarized noise';
        parameters.epoch_range           = [0.050 1.049];
        parameters.condition_names       = {   ...
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
        parameters.ITI                   = 11;
        
    case 3
        parameters.group                 = '1000ms binarized noised';
        parameters.epoch_range           = [0.050 1.049];
        parameters.condition_names               = {   ...
            'White Noise' ...
            'Binary White Noise' ...
            'Binary Pink Noise' ...
            'Binary Brown Noise' ...
            'Gratings(0.36 cpd)' ...
            'Gratings(0.73 cpd)' ...
            'Gratings(1.45 cpd)' ...
            'Gratings(2.90 cpd)' ...
            'Plaid'...
            'Blank'};
        parameters.ITI                   = 11;
        
    case 4
        parameters.group                 = '1000ms binarized noised';
        parameters.epoch_range           = [0.050 1.049];
        parameters.condition_names               = {   ...
            'House(High Contrast)' ...
            'House(Mid Contrast)' ...
            'House(Low Contrast)' ...
            'Face(High Contrast)' ...
            'Face(Mid Contrast)' ...
            'Face(Low Contrast)' ...
            'Pink Noise(High Contrast)' ...
            'Pink Noise(Mid Contrast)' ...
            'Pink Noise(Low Contrast)'...
            'Gratings(64cpi High Contrast)'...
            'Gratings(64cpi Mid Contrast)'...
            'Gratings(64cpi Low Contrast)'...
            'Blank'};
        parameters.ITI                   = 14;
end



parameters.baseline_condition = find(~cellfun(@isempty,strfind(condition_names, 'Blank')));

return