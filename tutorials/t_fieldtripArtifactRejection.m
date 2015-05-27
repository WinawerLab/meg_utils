% Tutorial from field trip on Visual artifact rejection
%
% http://www.fieldtriptoolbox.org/tutorial/visual_artifact_rejection
%
% Requires example data set, which can be found on Winawer Lab server or on
% fieldtrip's website:
%   /Volumes/server/Projects/MEG/code/ArtifactMEG.ds/PreprocData.mat
%   ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/ArtifactMEG.zip

% If the example data is available from Winawer Lab server, go there:
if exist('/Volumes/server/Projects/MEG/SampleData/PreprocData.mat', 'file')
    cd /Volumes/server/Projects/MEG/SampleData/
end


% Add fieldtrip path
field_trip_pth = '/Volumes/server/Projects/MEG/code/fieldtrip';
meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');

% Load the preprocessed data
load PreprocData dataFIC

%% To browse through the data trial by trial while viewing all channels write:
%  Click through the trials using the > button to inspect each trial.
%

cfg          = [];
cfg.method   = 'trial';
cfg.alim     = 1e-12; 
dummy        = ft_rejectvisual(cfg,dataFIC);

% Notice the weird vertical stripes. This is apparently due to a scale
% difference between MEG and EEG. To fix, see next cell:


%% If your dataset contains MEG and EEG channels (like this dataset), the
%  MEG and EEG channels are scaled differently when using only cfg.alim
%  (the EEG channels show up as big black bars on the screen). One of the
%  reasons to record EOG, EMG or ECG is to check these channels while
%  identifying eye, muscle and heart artifacts. 
%  
% The following code can be used to scale MEG and EEG channels both properly:

cfg          = [];
cfg.method   = 'trial';
cfg.alim     = 1e-12; 
cfg.megscale = 1;
cfg.eogscale = 5e-8;
dummy        = ft_rejectvisual(cfg,dataFIC);

% In trial 15 notice the slower drift observed over a larger group of
% sensors. This is most likely due to a head movement.

% Trial 84 shows an artifact which is caused by the electronics. Notice the
% jump in sensor MLT41:

% By browsing through the trials, related artifacts become evident (trial
% 15, 36, 39, 42, 43, 45 ,49, 50, 81, 82 and 84). They should be marked as
% 'bad'. After pressing the 'quit' button the trials marked 'bad' are now
% removed from the data structure.

