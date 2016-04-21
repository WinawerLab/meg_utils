function parameters = gamma_get_parameters(session_num)
%% parameters = gamma_get_experimental_parameters(session_num)
%  outputs a structure with the experimental parameters corresponding to
%  the given experiment

%% 
parameters = struct('sessionNumber',[],'epochRange',[],...
    'conditionVector', [],  'baselineCondition',[],'ITI',[],...
    'badChannels', [], 'badEpochs', [], 'envDenoise', [], ...
    'pcaDenoise', [], 'nBoot', [], 'f_sel', [] );

parameters.sessionNumber = session_num;
% becomes 1 if the data is denoised
parameters.envDenoise    = 0; 
parameters.pcaDenoise    = 0;


if session_num < 5
    condition_group = 1;
elseif ismember(session_num, 5:9)
    condition_group = 2;
elseif ismember(session_num, 10:16)
    condition_group = 3;
elseif session_num > 16 && session_num <99
    condition_group = 4;
end


switch condition_group
    case 1
        parameters.epochRange           = [0.050 0.549];
        parameters.ITI                   = [];
        parameters.baselineCondition    = 11;
        
    case 2
        parameters.epochRange           = [0.050 1.049];
        parameters.ITI                   = 11;
        parameters.baselineCondition    = 10;
        
    case 3
        parameters.epochRange           = [0.050 1.049];
        parameters.ITI                   = 11;
        parameters.baselineCondition    = 10;
        
    case 4
        parameters.epochRange           = [0.050 1.049];
        parameters.ITI                   = 14;
        parameters.baselineCondition    = 13;
end

return