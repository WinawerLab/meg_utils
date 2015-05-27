% Tutorial from field trip to remove EOG artifacts using ICA
%
% http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_eog_artifacts
%
% Requires example data set, which can be found on Winawer Lab server or on
% fieldtrip's website:
%   /Volumes/server/Projects/MEG/SampleData/ArtifactMEG.ds/
%   ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/ArtifactMEG.zip

% If the example data is available from Winawer Lab server, go there:
if exist('/Volumes/server/Projects/MEG/SampleData/ArtifactMEG.ds', 'dir')
    cd /Volumes/server/Projects/MEG/SampleData/
end


% Add fieldtrip path
field_trip_pth = '/Volumes/server/Projects/MEG/code/fieldtrip';
meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');


% preprocessing of example dataset
cfg = [];
cfg.dataset            = 'ArtifactMEG.ds';
cfg.trialdef.eventtype = 'trial';
cfg = ft_definetrial(cfg);

cfg.channel            = 'MEG';
cfg.continuous         = 'yes';
data = ft_preprocessing(cfg);

% downsample the data to speed up the next step
cfg = [];
cfg.resamplefs = 300;
cfg.detrend    = 'no';
data = ft_resampledata(cfg, data);

% perform the independent component analysis (i.e., decompose the data)
cfg        = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB

comp = ft_componentanalysis(cfg, data);


% Redefine cfg
cfg = [];
cfg.method = 'wavelet';
    
% Get either the Fourier spectrum 
freq = ft_freqanalysis(cfg, comp);

% ... or timelocked averaged events
cfg = [];
timelock = ft_timelockanalysis(cfg, comp);

% plot the components for visual inspection
figure
cfg = [];
cfg.component = [1:20];       % specify the component(s) that should be plotted
cfg.layout    = 'CTF151.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

cfg = [];
cfg.layout = 'CTF151.lay'; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
ft_databrowser(cfg, comp)

% remove the bad components and backproject the data
cfg = [];
cfg.component = [9 10 14 24]; % to be removed component(s)
data = ft_rejectcomponent(cfg, comp, data);