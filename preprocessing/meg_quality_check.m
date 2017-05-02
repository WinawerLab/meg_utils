%% MEG Quality check script

% This is a script to check some basics of, and possible artifacts in your MEG
% data. 

%% 0. Set paths and other settings

% Some default settings
verbose = true;
samplesToUse = 200000;%[]; % How many samples do you want to use for the quality check, leave empty if you want whole data set

% Path to dataset
pth = '/Volumes/server/Projects/MEG/Gamma/Data/';
session = '20_Gamma_3_22_2016_subj010';
sqdfile = 'R0852_Gamma_3.22.16.sqd';

dataPth = fullfile(pth,session,'raw',sqdfile);

% Add fieldtrip to path with the ToolboxToolbox
if ~exist('ft_read_header','file')
    tbUse('Fieldtrip');
end

% In case of Yokogawa dataset
refChannels = 158:160;
eyeChannels = 177:179;
dataChannels = 1:157;
triggerChannels = 161:168;
photoDiodeChannel = 192;

%% 1. Load data and header information

% Load header
hdr = ft_read_header(dataPth);

% Load data
data = ft_read_data(dataPth);


%% 2. Get some basics
results = [];
results.fs        = hdr.Fs;
results.nsamples  = hdr.nSamples;
results.units     = hdr.chanunit{1};
results.nchannels = size(data,1);
results.datapth   = dataPth;

if verbose
    fprintf('\n')
    fprintf('Some basics.. \n')
    fprintf('Recording frequency : %d Hz\n', results.fs)        % Recording frequency
    fprintf('Nr of samples       : %d\n', results.nsamples)           % Recording duration
    fprintf('Data are in units of: %s\n',results.units) % Units of data
    fprintf('Total nr of channels: %d\n', results.nchannels)
end

% Take a sample, or not
if ~isempty(samplesToUse>1)
    dataAll = data;
    data = data(:,1:samplesToUse);
    results.sampleSize = samplesToUse;
end

%% 3. Plot timeseries of data

if verbose
    % Show sample timeseries channels
    figure;
    subplot(5,1,1); hold all; title('Data Channels')
    for ii = dataChannels; plot(data(ii,:)); end
    
    subplot(5,1,2); hold all; title('Environmental reference channels')
    for ii = refChannels; plot(data(ii,:)); end
    
    subplot(5,1,3);hold all; title('Eyetracking channels')
    for ii = eyeChannels; plot(data(ii,:)); end
    
    subplot(5,1,4); hold all; title('Trigger channels')
    for ii = triggerChannels; plot(data(ii,:)); end
    
    subplot(5,1,5); hold all; title('Photodiode channel')
    for ii = photoDiodeChannel; plot(data(ii,:)); end
end

% Triggers send?
trigger    = meg_fix_triggers(dataAll(triggerChannels,:)');
results.triggerVals  =  unique(trigger(trigger~=0));
results.triggerTotal = sum(find(trigger(trigger~=0)));

for trig = results.triggerVals'
    results.triggerCount(trig) = length(find(trigger==trig));
end

if verbose
    figure; histogram(trigger(trigger~=0)); xlabel('trigger values'); ylabel('frequency')
    fprintf('Following triggers found: %d\n', results.triggerVals)
    fprintf('Amount per trigger: %d \n', results.triggerCount)
    fprintf('Total amount of triggers found: %d\n', results.triggerTotal)
end


%% 4. Plot spectrum
numFreq     = size(data, 2);
t           = (1:numFreq)/results.fs;
f           = (0:length(t)-1)/max(t);
F           = (abs(fft(data,[],2))/length(t)*2);

[pks, locs] = findpeaks(smooth(F(1,:),20),f,'MinPeakProminence', 3E-15,'MinPeakDistance',1);
locs = locs(locs>0.1); % eliminate frequencies close to 0
if verbose
    figure; 
    cla; subplot(211)
    plot(f,F(1,:)); hold on;
    plot(f(round(locs)),F(1,round(locs)),'ro');
    title('Spectrum of data samples from channel 1')
    set(gca,'YScale','log','XLim',[0 250]);
    xlabel('Frequency [Hz]'); ylabel('Amplitude (pT)')
    
    subplot(212)
    plot(f,F(25,:));
    title('Spectrum of data samples from channel 25')
    set(gca,'YScale','log','XLim',[0 250]);
    xlabel('Frequency [Hz]'); ylabel('Amplitude (pT)')
    
end



%   Filters (e.g., high-pass, line noise) -- Power spectrum density (Welch), average/smooth??     find peaks?? how to determine
%   a high-pass filter

% Power spectra of 3 reference channels (158:160)
%   Best fit 1/f?
%   Line noise
%   Other large peaks?
%   Jumps???
%
% How much line noise?
%   the power (relative to what??)
%   the distribution across channels and time
%   which harmonics?
%
% Jumps, artifacts, etc
% cfg.dataset = dataPth;
% cfg.continuous = 'yes';
% cfg.trl = [];

cfg = [];

conditions = trigger(find(trigger));
cfg.dataset = dataPth;
cfg.trialdef.trig = dataAll(triggerChannels,:);
cfg.trialdef.eventtype = 'trial';
cfg.trialdef.pre = 0;
cfg.trialdef.post = 1;
cfg.trialdef.trialFunHandle = @ssmeg_ft_trial_fun;
cfg.trialdef.conditions = conditions;
cfg.trialdef.onsets = trigger;

[trl, Events] = ssmeg_ft_trial_fun(cfg, trigger, conditions);

cfg = [];
cfg.trl = trl;
cfg.continuous = 'no';

[cfg, artifact] = ft_artifact_jump(cfg,dataPth);


% Other weird peaks in data
%
% Saturated channels
%   which, when, how many?
%
% Save out images of power spectrum and spectrogram of every channel
%% Look at empty room data
pth = '/Volumes/server/Projects/MEG/Gamma_BR/emptyRoomData/Empty_Room_Gamma Settings_6.2.2016_1135am.sqd';
hdr = ft_read_header(pth);
data = ft_read_data(pth);



% also consider: http://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsFilter