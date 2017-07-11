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
if ~exist('dataPth','var') || isempty(dataPth)
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

%% Load data and header information

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

%% Get some basics and write them to log file

results = [];
results.fs        = hdr.Fs;
results.nsamples  = hdr.nSamples;
results.units     = hdr.chanunit{1};
results.nchannels = size(data,1);
results.datapth   = dataPth;

logfile = fullfile(resultsDir, 'logfile.txt');
fid = fopen(logfile, 'w');

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


%% Derived parameters

% convert units to picotesla
switch lower(results.units)
    case 't',        toPT = 1E12;
    case {'nt' 'n'}, toPT = 1E3;
    case {'pt' 'p'}, toPT = 1;
    case {'ft' 'f'}, toPT = 1E-3;
    otherwise, error('Unrecognized unit %s', results.units);
end
       

% Select the range of the dataset for ananlysis
if isempty(samplesToUse), samplesToUse = size(data, 2); end
samplesToUse = 1:min(samplesToUse, size(data,2));

sampleData = data(:,samplesToUse);
results.sampleSize = length(samplesToUse);

%% Plot sample timeseries of data
t = (1:length(samplesToUse)) / results.fs;

% Show sample timeseries channels
% fH =  megNewGraphWin([],'tall'); 
fH = figure;
set(fH, 'Name', 'Sample Time Series'); 

subplot(5,1,1); hold all; title('Data Channels')
plot(t, sampleData(dataChannels,:)*toPT);
ylabel('Magnetic Field (pT)');

subplot(5,1,2); hold all; title('Environmental reference channels')
plot(t, sampleData(refChannels,:)*toPT);
ylabel('Magnetic Field (pT)');

subplot(5,1,3);hold all; title('Eyetracking channels')
plot(t, sampleData(eyeChannels,:));
ylabel('Position (mm??)');

subplot(5,1,4); hold all; title('Trigger channels')
plot(t, sampleData(triggerChannels,:));

subplot(5,1,5); hold all; title('Photodiode channel')
plot(t, sampleData(photoDiodeChannel,:));
xlabel('Time (seconds)');

hgexport(fH, fullfile(resultsDir, 'SampleTimeSeries.eps'));

close(fH)

%% Triggers 

% We use several trigger channels, each of which is binary, and the binary
% values are then converted to base 10. Before converting to base 10, the
% function checks for slight timing errors, for example if two trigger
% channels have a signal one sample apart, they are assumed to be
% synchronous. 

% the variable trigger is a time series of trigger values, equal in length
% to the MEG data set
trigger    = meg_fix_triggers(data(triggerChannels,:)');

% unique base 10 trigger values
results.trigger.vals  =  unique(trigger(trigger~=0));

for trig = results.trigger.vals'
    results.trigger.count(trig==results.trigger.vals) = sum(trigger==trig);
end

% total number of triggers detected 
results.trigger.total = sum(results.trigger.count);


% fH = megNewGraphWin([], [], 'Name', 'Triggers'); 
fH = figure;
histogram(trigger(trigger~=0)); xlabel('trigger values'); ylabel('frequency')
hgexport(fH, fullfile(resultsDir, 'TriggerHistogram.eps'));
close(fH);

fprintf(fid, '\n\nTRIGGERS\n');
fprintf(fid, 'TriggerVal\tCount\n');
for ii = 1:length(results.trigger.vals)
    fprintf(fid, '%d\t%d\n', results.trigger.vals(ii), results.trigger.count(ii));
end
fprintf(fid, 'Total amount of triggers found: %d\n', results.trigger.total);



%% CONTINUE HERE!!!
% CONTINUE HERE!!!
%% CONTINUE HERE!!!

%% 4. Plot spectrum
numFreq     = size(data, 2);
t           = (1:numFreq)/results.fs;
f           = (0:length(t)-1)/max(t);
F           = abs(fft(data,[],2))/length(t)*2;

[pks1, locs1] = findpeaks(smooth(F(1,:),20),f,'MinPeakProminence', 8E-16,'MinPeakDistance',1);
locs1 = locs1(locs1>1); % eliminate frequencies close to 0
pks1 = pks1(locs1>1);

fH = figure;
cla; subplot(211)
plot(f,F(1,:)); hold on;
plot(locs1,pks1,'ro');
title('Spectrum of data samples from channel 1')
set(gca,'YScale','log','XLim',[0 250]);
xlabel('Frequency [Hz]'); ylabel('Amplitude (pT)')

[pks25, locs25] = findpeaks(smooth(F(1,:),20),f,'MinPeakProminence', 8E-16,'MinPeakDistance',1);
locs25 = locs25(locs25>1); % eliminate frequencies close to 0
pks25 = pks25(locs25>1);

subplot(212)
plot(f,F(25,:)); hold on;
plot(locs25,pks25,'ro');
title('Spectrum of data samples from channel 25')
set(gca,'YScale','log','XLim',[0 250]);
xlabel('Frequency [Hz]'); ylabel('Amplitude (pT)')
hgexport(fH, fullfile(resultsDir, 'SpectrumChannel1_25.eps'));

results.peak_locs = unique([locs1, locs25]);

fprintf(fid, 'Following sharp peaks found in spectrum, in Hz: %d\n', results.peak_locs);

% To DO
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


% Other weird peaks in data
%
% Saturated channels
%   which, when, how many?
%
% Save out images of power spectrum and spectrogram of every channel
%
%% Define 1 second trials in FT struct in order to get SQD jumps
cfg = [];

cfg.dataset             = fullfile(dataPth);
cfg.trialdef.trig       = triggerChannels;
cfg.trialdef.eventtype  = 'trial';
cfg.trialdef.pre        = 0;
cfg.trialdef.post       = 1;
cfg.trialdef.trialFunHandle = @ssmeg_ft_trial_fun;

[conditions, onsets] = find(trigger);

cfg.trialdef.conditions    = conditions;
cfg.trialdef.onsets       = onsets;

[trl, Events] = ssmeg_ft_trial_fun(cfg, onsets, conditions);
cfg.trl       = trl;
cfg.event     = Events;

%% Get trials with SQD jumps
cfg                     = [];
cfg.trl                 = trl;
cfg.continuous          = 'yes';
cfg.datafile            = dataPth;
cfg.headerfile          = dataPth;
data_hdr                = load('hdr'); data_hdr = data_hdr.hdr;
cfg.layout              = ft_prepare_layout(cfg, data_hdr);

% Check for MEG channels, define cutoff and add no padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 20;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';

% make the process interactive or not
cfg.artfctdef.zvalue.interactive = 'no';

% Takes time
[cfg, artifact_jump] = ft_artifact_zvalue(cfg);

% Report artifact jumps
fprintf(fid, 'Following artifact (squid) jumps found using 1 sec epochs: %d\n, using a zvalue threshold of %d', artifact_jump, cfg.artfctdef.zvalue.cutoff);



fclose(fid);






%% Look at empty room data
dataPth = '/Volumes/server/Projects/MEG/Gamma_BR/emptyRoomData/Empty_Room_Gamma Settings_6.2.2016_1135am.sqd';
hdr = ft_read_header(dataPth);
data = ft_read_data(dataPth);



% also consider: http://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsFilter