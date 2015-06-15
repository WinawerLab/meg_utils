%% SSEEG HPC script to denoise
function HPC_SSEEG(subjects)


%% Get variables

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = false;


%% Set analysis variables
rootPath                      = which('HPC_SSEEG');
rootPath                      = fileparts(rootPath);

project_pth                   = fullfile(rootPath,'Data');

%% Prepare and solve GLM

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
epoch_time = [0 0.9940];
s_rate_eeg = 1000;
T = diff(epoch_time)+ 1/s_rate_eeg;
freq = megGetSLandABfrequencies((0:150)/T, T, 12/T);



% denoise parameters (see denoisedata.m)
opt.pcchoose          = 1.05;
opt.npcs2try          = 10;
opt.npoolmethod       = {'r2','n',40};
opt.verbose           = true;
opt.pcn               = 10;
opt.savepcs           = 0;
optsl = opt;
optbb = opt;
optbb.preprocessfun   = @(x)hpf(x, freq.sl);       % preprocess data with a high pass filter for broadband analysis
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband
npc_used              = [20,25,30,35,40,45,50,60];
allResults            = [];
allEvalout            = [];

% The denoise algorithm needs:
% data      : time series [channel x time samples x epoch]
% design    : design matrix [epoch x nconds]
for whichSubject = subjects
    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(project_pth, 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(project_pth, 's0%d_design.mat'),whichSubject)); design = tmp.design;
    
    % Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    %% Denoise the broadband data
    for nr_pc = 1:length(npc_used)
        optbb.pcstop = -npc_used(nr_pc);
        fprintf('\tnpcs = %d\n', npcs(nr_pc));
        
        [results, evalout] = denoisedata(design,sensorData,noisepooldef,evalfun,optbb);
        allResults{nr_pc} = results;
        allEvalout{nr_pc} = evalout;
        
        clear results; clear evalout;
    end
    
    
    fname = sprintf(fullfile(project_path, 's0%d_denoisedData_varynpcs'),whichSubject);
    parsave([fname '_bb.mat'], 'allResults', allResults, 'allEvalout', allEvalout, 'opt', optbb)
    
    
    %%  Denoise for stimulus-locked analysis
    allResults            = [];
    allEvalout            = [];
    
    
    for nr_pc = 1:length(npc_used)
        optbb.pcstop = -npc_used(nr_pc);
        fprintf('\tnpcs = %d\n', npcs(nr_pc));
        
        [results, evalout] = denoisedata(design,sensorData,noisepooldef,evalfun,optsl);
        allResults{nr_pc} = results;
        allEvalout{nr_pc} = evalout;
        
        clear results;
    end
    fname = sprintf(fullfile(project_path, 's0%d_denoisedData_varynpcs'),whichSubject);
    parsave([fname '_sl.mat'], 'allResults', allResults, 'allEvalout', allEvalout, 'opt', optsl)
    
end