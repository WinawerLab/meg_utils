%% s_GAMMA_data_pipeline
%  new wrapper script for a new analysis pipeline
%  goals:      ~ eliminate redundancies (processes that are repeated
%              in both the analysis and denoising scripts)
%              ~ allow both the analysis and visualization functions to
%              take either denoised or undenoised data as input
%              ~ save files in organized structures and locations
%              ~ allow for parallel processing (and eventual HPC
%              compatability
%
% Nicholas Chua

%% Paths and options

% add field trip paths
if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end


whichSessions          = 18;
environmentalChannels  = 158:160;
dataChannels           = 1:157;

environmentalDenoising = false;
PCADenoising           = true;
nBoot                  = 5;

%% Loop over sessions
for sessionNum = whichSessions
    %% get paths/parameters
    
    sessionPath = meg_gamma_get_path(sessionNum);
    params = gamma_get_parameters(sessionNum);
    
    params.nBoot =  nBoot;
    
    %    dataPath    = fullfile(sessionPath, 'raw');
    %    savePath    = fullfile(path_to_data, 'processed');
    %
    %    if ~exist(savePath, 'dir'), mkdir(savePath); end
    %% preprocessing
    
    % badChannels and badEpochs to be preserved before env denoise
    [tmp, conditionVector, params.badChannels,...
        params.badEpochs] = gamma_preprocess_raw(sessionNum);
    
    %% Denoising
    
    % denoise using environmental channels
    if environmentalDenoising
        tmp               = meg_environmental_denoising(tmp,...
                             environmentalChannels, dataChannels);
        params.envDenoise = 1;
    end
    
    % remove badEpochs and badChannels from ts after environmental
    % denoising
    ts                            = tmp(:, ~params.badEpochs, :); 
    clear tmp;
    ts(:,:,params.badChannels)    = NaN;
    params.conditionVector        = conditionVector(~params.badEpochs);
    
    % prepare timeseries for PCA and spectral analysis by removing
    % ITI epochs from the timeseries and conditionVector
    
    ITIepochs  = params.conditionVector == params.ITI;
    ts         = ts(:, ~ITIepochs, :);
    params.conditionVector = params.conditionVector(params.conditionVector ~= params.ITI);
    
    % denoise using Kuper's PCA
    if PCADenoising
        ts                = gamma_denoise(ts, params);
        params.pcaDenoise = 1;
    end
    %% Spectral analysis, bootstraping, and fitting
    
    % results is a struct with fields:
    % spectral_data_mean, fit_bl_mn, w_pwr_mn, w_gauss_mn, ...
    % gauss_f_mn, fit_f2_mn
    results = gamma_spectral_analysis(ts(:,:,dataChannels), params);
    
    
end



