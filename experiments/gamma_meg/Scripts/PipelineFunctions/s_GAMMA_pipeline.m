%% s_GAMMA_data_pipeline
%  new wrapper script for a new analysis pipeline
%  principles: ~ eliminate operation redundency (processes that are repeated
%              in both the analysis and denoising scripts
%              ~ allow both the analysis and visualization functions to
%              take denoised and undenoised data
%              ~ save files in an organized structure
%              ~ allow for parallel processing (and eventual HPC
%              compadibility
%
% Nicholas Chua

%% Paths and options

% add field trip paths
if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end

whichSessions = 18;
denoiseData   = true;
environmentalChannels = 158:160;
dataChannels = 1:157;

environmentalDenoising = true;

%% Loop over sessions
for sessionNum = whichSessions
   %% get paths to data and save locations
   
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
   ts = tmp(:, ~badEpochs, ~badChannels);
   conditionVector = conditionVector(~badEpochs);
   
   % denoise using Kuper's PCA
   denoisedTs = gamma_denoise(ts, conditionVector, sessionNum);
   
   %% Spectral analysis, bootstraping, and fitting
   
    
    
end

