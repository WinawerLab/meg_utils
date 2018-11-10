%% t_forwardmodels

% This is a tutorial to describe general steps to use MEG forward models
% and predict responses from vertices on the surface to sensors
%          
%      https://osf.io/rtseb/?action=download

% Completed before running this tutorial:
%   0.1 Download data from OSF, move to data folder and unzip
%   0.2 Run FreeSurfer's recon-all on subject T1 anatomy
%   0.3 Set up new Brainstorm protocol, add subject anatomy folder and data
%       (Fiducials can be found in downloaded folder)

% Overview:
%   1. Set up paths, project and subject
%   2. Run benson docker if V1-V3 templates don't exist
%   3. Downsample V1-V3 templates  to Brainstorm mesh
%   4. Load gain matrix 
%   5. Get V1 template
%   6. Create a prediction from V1 vertices to MEG sensors
%   7. Plot the mean predicted fourier amplitudes and brainstorm mesh


% Dependencies:
%   * Forward model synchrony project repository (https://github.com/elinekupers/fmsproject)
%   * Brainstorm (https://github.com/brainstorm-3/)
%   ---> NB: Brainstorm GUI has to be open

% see also: t_preprocessSampleData.m


%% 1. Set up paths, project and subject

% Define Brainstorm and FreeSurfer data base folder
bsDB            = '/Volumes/server/Projects/MEG/brainstorm_db/';
fsDB            = '/Volumes/server/Freesurfer_subjects/';


% Protocol name in Brainstorm database
protocolName     = 'MEGSampleData';

% Subject name in Brainstorm project
subjectBS         = 'wlsubj004';

% Subject name in FreeSurfer 
subjectFS         = 'wlsubj004';


% Get data and anatomy folder in brainstorm database
d = dir(fullfile(bsDB, protocolName, 'data', subjectBS, 'R*'));

% If the data folder can't be found, spit out and error:
if size(d,1) == 0
    error('Sample dataset does not exist in brainstorm yet. Open GUI by typing ''brainstorm'' in MATLAB command window and add protocol with data')
end

dataDir = fullfile(d(1).folder, d(1).name);    
anatDir = fullfile(bsDB, protocolName, 'anat', subjectBS);

    
%% 2. Run benson docker if V1-V3 templates don't exist

fsFolder = fullfile(fsDB, subjectFS);
fsDirAngle = dir(fullfile(fsFolder, 'surf', '*template*angle.mgz'));
fsDirEccen = dir(fullfile(fsFolder, 'surf', '*template*eccen.mgz'));
fsDirArea = dir(fullfile(fsFolder, 'surf', '*template*area.mgz'));


if size(fsDirAngle,1) == 0
    str = sprintf('docker run -ti --rm -v %s:/input \\nben/occipital_atlas:latest', fsFolder);
    system(str)
    
    fsDirAngle = dir(fullfile(fsFolder, 'surf', '*template*angle.mgz'));
    fsDirEccen = dir(fullfile(fsFolder, 'surf', '*template*eccen.mgz'));
    fsDirArea  = dir(fullfile(fsFolder, 'surf', '*template*area.mgz'));
end

%% 3. Downsample V1-V3 templates  to Brainstorm mesh

bsDirAngle = dir(fullfile(bsDB, protocolName, subjectBS, 'anat', '*template*angle.mat'));

if size(bsDirAngle,1) == 0
    % Add matlab compatible freesurfer code
    if ~exist('MRIRead')
        addpath(genpath('/Applications/freesurfer/matlab'));
        addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
    end
    
    % Get templates from FS, and create downsampled BS templates
    interp_retinotopy(bsDB, fsDB, subjectFS, subjectBS, protocolName)
end
    

%% 4. Load gain matrix 

% Define vector that can truncate number of sensors 
G_constrained = getGainMatrix(dataDir);

%% 5. Get V1 template

% Get V1 template limited to 11 degrees eccentricity. 
template = getTemplate(anatDir, 'V1', 11);

%% 6. Create a prediction from V1 vertices to MEG sensors

% Simulate coherent and incoherent source time series and compute
% predictions from forward model (w)
nrTimePoints   = 1000; % ms (1000 Hz sample rate)
nrEpochs       = 1;    % 
freq           = 1;    % frequency of simulated vertex sine wave
predictions    = getForwardModelPredictions(G_constrained, template.V1StimEccen, freq, nrTimePoints, nrEpochs);

%% 7. Plot the mean predicted fourier amplitudes and brainstorm mesh

% Take get fourier amplitudes
amps.c = abs(fft(predictions.c,[],2));
amps.i = abs(fft(predictions.i,[],2));

% Take the mean across epochs
w.V1c = mean(amps.c(:,2,:),3);
w.V1i = mean(amps.i(:,2,:),3);

% Visualize brainstorm mesh
colors = jet(3);
visualizeBrainstormMesh(anatDir, colors)

% Visualize predictions from forward model
figure(98); clf;
subplot(121)
megPlotMap(w.V1c,[],[],'bipolar')
subplot(122)
megPlotMap(w.V1i,[],[],'bipolar')
