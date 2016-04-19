%% s_GAMMA_data_pipeline
%  new wrapper script for a new analysis pipeline
%  goals:      ~ eliminate redundancies (processes that are repeated
%              in both the analysis and denoising scripts)
%              ~ allow both the analysis and visualization functions to
%              take denoised and undenoised data
%              ~ save files in organized structures
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
PCADenoising           = false;
nBoot                  = 5;

%% Loop over sessions
for sessionNum = whichSessions
    %% get paths/file suffix
    
    sessionPath = meg_gamma_get_path(sessionNum);
    
    %    dataPath    = fullfile(sessionPath, 'raw');
    %    savePath    = fullfile(path_to_data, 'processed');
    %
    %    if ~exist(savePath, 'dir'), mkdir(savePath); end
    %% preprocessing
    
    % badChannels and badEpochs to be preserved before env denoise
    [tmp, conditionVector, badChannels, badEpochs] = gamma_preprocess_raw(sessionNum);
    
    %% Denoising
    
    % denoise using environmental channels
    if environmentalDenoising
        tmp = meg_environmental_denoising(tmp, environmentalChannels, dataChannels);
    end
    
    % remove badEpochs and badChannels from ts after environmental
    % denoising
    ts = tmp(:, ~badEpochs, :); clear tmp;
    ts(:,:,badChannels) = NaN;
    conditionVector = conditionVector(~badEpochs);
    
    % denoise using Kuper's PCA
    if PCADenoising
        ts = gamma_denoise(ts, conditionVector, sessionNum);
    end
    %% Spectral analysis, bootstraping, and fitting
    
    % results is a struct with fields:
    % spectral_data_mean, fit_bl_mn, w_pwr_mn, w_gauss_mn, ...
    % gauss_f_mn, fit_f2_mn
    results = gamma_spectral_analysis(ts, conditionVector, nBoot, sessionNum);
    
    
end



