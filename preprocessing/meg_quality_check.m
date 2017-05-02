function meg_quality_check(dataPth, varargin)
%Function to check MEG data quality
% meg_quality_check(dataPth, varargin)
% 
% This is a script to check some basics of, and possible artifacts in your
% MEG data.
%
% Example 1: NYU Yokogawa 160
%  dataPth = '/Volumes/server/Projects/MEG/Gamma/Data/20_Gamma_3_22_2016_subj010/raw/R0852_Gamma_3.22.16.sqd';
%  meg_quality_check(dataPth, 'verbose', true, 'samplesToUse', 200000);
%
% Example 2: NYU-AD Yokogawa 224
%  dataPth = '/Volumes/server/Projects/MEG/NYUAD/0167_SSMEG_20170123/raw/0167_SSMEG_20170123_01.con';
%  meg_quality_check(dataPth, 'verbose', true, 'samplesToUse', []);



% Check inputs
if ~exist(dataPth,'var') || isempty(dataPth)
    error('Path needed for MEG data');
end

if ~exist(dataPth, 'file')
    error('Path %s does not exist', dataPth);
end

% optional inputs
for ii = 1:2:length(varargin)
   param = varargin{ii};
   val   = varargin{ii+1};
   switch lower(param)
       case 'verbose' 
           verbose = val;
       case 'samplestouse'
           samplesToUse = val;
   end
    
end

% Check for field trip
if ~exist('ft_read_header','file')
    error('Field trip is needed on your path.')
end

% Set some defaults
if ~exist('verbose', 'var'), verbose = true; end
if ~exist('samplesToUse', 'var'), samplesToUse = []; end

% results directory
resultsDir = fullfile(fileparts(dataPth), 'dataQuality');
mkdir(resultsDir);

%% 1. Load data and header information

% Load header
hdr = ft_read_header(dataPth);

% Load data
data = ft_read_data(dataPth);


% Check where data is from
dataFormat = ft_filetype(dataPth);
    
switch dataFormat
    case 'yokogawa_ave'
        % In case of Yokogawa dataset
        dataChannels        = 1:157;
        refChannels         = 158:160;
        eyeChannels         = 177:179;
        triggerChannels     = 161:168;
        photoDiodeChannel   = 192;
        
    case 'yokogawa_con'
        dataChannels        = 1:208;
        refChannels         = 209:211;
        triggerChannels     = 224:231;
        eyeChannels         = [];
        photoDiodeChannel   = 233;

end

%% 2. Get some basics
results = [];
results.fs        = hdr.Fs;
results.nsamples  = hdr.nSamples;
results.units     = hdr.chanunit{1};
results.nchannels = size(data,1);
results.datapth   = dataPth;

logfile = fullfile(resultsDir, 'logfile.txt');
fid = fopen(logfile, 'w');

if verbose
    fprintf(fid, '\n');
    fprintf(fid, 'Some basics.. \n');
    fprintf(fid, 'Recording frequency:\t%d Hz\n', results.fs);        % Recording frequency
    fprintf(fid, 'Nr of samples:\t%d\n', results.nsamples);           % Recording duration
    fprintf(fid, 'Data are in units of:\t%s\n',results.units); % Units of data
    fprintf(fid, '\n\n-----Channel Information----\n');
    fprintf(fid, 'Total nr of channels:\t%d\n', results.nchannels);
    fprintf(fid, 'Data channels:\t%d:%d\n',          min(dataChannels), max(dataChannels));
    fprintf(fid, 'Reference channels:\t%d:%d\n',     min(refChannels), max(refChannels));
    fprintf(fid, 'Trigger channels:\t%d:%d\n',       min(triggerChannels), max(triggerChannels));
    fprintf(fid, 'Eye channels:\t%d:%d\n',           min(eyeChannels), max(eyeChannels));
    fprintf(fid, 'Photodiode channel:\t%d\n',        photoDiodeChannel);

    
end
fclose(fid);

% Take a sample, or not
if isempty(samplesToUse), samplesToUse = size(data, 2); end
samplesToUse = min(samplesToUse, size(data,2));

sampleData = data(:,1:samplesToUse);
results.sampleSize = samplesToUse;


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
dataPth = '/Volumes/server/Projects/MEG/Gamma_BR/emptyRoomData/Empty_Room_Gamma Settings_6.2.2016_1135am.sqd';
hdr = ft_read_header(dataPth);
data = ft_read_data(dataPth);



% also consider: http://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsFilter