function s_GAMMA_pipeline(whichSessions)
%  
% This is a script to run separate parts of analysis pipeline of the MEG 
% Gamma project. One can save results or timeseries of every part, so that
% this part of the analysis can be skipped in the future.
%
% This pipeline workflow:
%   (0) Load toolboxes and define options
%   (1) Loading SQD files and preprocessing timeseries
%   (2) Denoise timeseries (environmentally or MEG Denoise)
%   (3) Compute spectra and fit Broadband / Gamma model
%   (4) Visualize data

% Session descriptions (more details in meg_get_parameters.m):
%       Session 2:9         500 ms stimuli, no the binary pink and brown noise stimuli
%       Session 5:9         1000 ms stimuli, can be denoised
%       Session 10:16       1000 ms stimuli, binary noise stimuli, can be denoised
%       Session 13          noisy and should not be used
%       Session 17-22       faces, houses, gratings, noise patterns at low, medium, and high contrast

% Dependencies:
%       MEG utils toolbox   (for loading, preprocessing, visualing data)
%       Fieldtrip toolbox   (for loading, visualing data)
%       MEG denoise toolbox (for preprocessing and denoising data)
%
% First version: Nicholas Chua (April 2016)
%       7.5.2016: Clean up (EK)
%       7.7.2016: Adding option to run script on HPC (EK)

%% ------------------------------------------------------------------------
%  -------------------- 0. Define paths and options -----------------------
%  ------------------------------------------------------------------------

% Do you want to run the analysis on the HPC?
opt.HPC = true;

if opt.HPC
    % Define root path
    rootPath = which(mfilename);
    rootPath = fileparts(fileparts(fileparts(rootPath)));   
    % Add denoiseproject
    addpath(genpath('/scratch/ek99/denoiseproject'));    
else % Add toolboxes located on the desktop
    if isempty(which('ft_prepare_layout'))
        meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
    end

    if isempty(which('dfdPreprocessData'))
        run '~/matlab/git/denoiseproject/dfdAddPaths'
    end

    if isempty(which('meg_gamma_get_path'))
        addpath(genpath('~/matlab/git/meg_utils'))
    end   
end

% Which sessions do you want to analyze?
% whichSessions          = 19;

% Parameters for preprocessing 
opt.dataChannels           = 1:157;
opt.environmentalChannels  = 158:160;
opt.triggerChannels        = 161:164;
opt.environmentalDenoising = true;          % use 3 environm. channels to denoise?
opt.MEGDenoise             = true;         % use MEG Denoise?
opt.nBoot                  = 100;             % nr of bootstraps
opt.fs                     = 1000;          % sampling rate (milliseconds)
% Epoch start and end is hardcoded in gamma_get_parameters(sessionNum) <-- do we want this?
% opt.epochStartEnd          = [0.050 1.049]; % start and end of epoch (seconds) 
opt.varThreshold           = [0.05 20];     % Threshold for variance in order to define noisy data
opt.badChannelThreshold    = 0.2;           % Threshold for proportion of epochs to be bad before eliminating whole channel
opt.badEpochThreshold      = 0.2;           % Threshold for proportion of channels to be bad before eliminating whole epoch

% General parameters for saving and printing
opt.verbose                = true;
opt.saveData               = true;
opt.saveFigures            = false;


%% Loop over sessions
for sessionNum = whichSessions
    
    
    %% --------------------------------------------------------------------
    %  --------------------- 1. Load data and preprocess ------------------
    %  --------------------------------------------------------------------
    
    % Get session paths and parameters 
    if opt.HPC
        projectPath = fullfile(rootPath,'HPC','Data');
        dataPth    = sprintf('%02d_Gamma_*subj*', sessionNum);
        d = dir(fullfile(projectPath, dataPth));
        sessionPath = fullfile(projectPath, d.name);   
        disp(sessionPath)
    else sessionPath = meg_gamma_get_path(sessionNum); 
    end
    
    params = gamma_get_parameters(sessionNum);
    opt.params      = params; clear params;
    opt.sessionPath = sessionPath; clear sessionPath;
    
    % Preprocessing (badChannels and badEpochs to be preserved before env
    % denoise)
    [ts, opt] = gamma_preprocess_raw(sessionNum, opt);
    
    %% --------------------------------------------------------------------
    %  ----------------------- 2. Denoise timeseries ----------------------
    %  --------------------------------------------------------------------
    
    % Remove badEpochs and badChannels from ts after environmental denoising
    ts                              = ts(:, ~opt.params.badEpochs, opt.dataChannels); 
    ts(:,:,opt.params.badChannels)  = NaN;
    opt.params.conditions           = opt.params.conditions(~opt.params.badEpochs);
    
    % prepare timeseries for PCA and spectral analysis by removing
    % ITI epochs from the timeseries and conditionVector 
    ITIepochs                       = opt.params.conditions == opt.params.ITI;
    ts                              = ts(:, ~ITIepochs, :);
    opt.params.conditions           = opt.params.conditions(opt.params.conditions ~= opt.params.ITI);
    
    % Denoise using MEG Denoise is requested
    if opt.MEGDenoise
        ts = ts(:,:,~opt.params.badChannels);        
        ts = gamma_denoise(ts, opt);

        sz = size(ts);
        ts = reshape(ts, [], sz(3));
        ts = to157chan(ts, ~opt.params.badChannels, 'nans');
        ts = reshape(ts, sz(1), sz(2), []);     
    end
    
    %% --------------------------------------------------------------------
    %  -------------------- 3. Spectral analyses --------------------------
    %  --------------------------------------------------------------------
    
    % Spectral analysis, bootstraping, and fitting
    [results, opt] = gamma_spectral_analysis(ts(:,:,opt.dataChannels), opt);
    
    %% --------------------------------------------------------------------
    %  ---------------------- 4. Visualization ----------------------------
    %  --------------------------------------------------------------------
    if ~opt.HPC; gamma_visualize_spectral(results, opt); end
    
end



