% Tutorial from field trip on artifact rejection
%
% http://www.fieldtriptoolbox.org/example/use_independent_component_analysis_ica_to_remove_eog_artifacts
%
% Requires example data set, which can be found on Winawer Lab server or on
% fieldtrip's website:
%   /Volumes/server/Projects/MEG/code/ArtifactMEG.ds/
%   ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/ArtifactMEG.zip

% If the example data is available from Winawer Lab server, go there:
if exist('/Volumes/server/Projects/MEG/code/', 'dir')
    cd /Volumes/server/Projects/MEG/code/
end

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

TO BE CONTINUED