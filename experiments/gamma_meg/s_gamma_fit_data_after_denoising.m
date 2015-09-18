%% Script to model gamma and broadband at the same time, after denoising

% Define variables
subjects           = 9;
fs                = 1000;
nboot             = 10;
trigger_channels  = 161:164;
data_channels     = 1:157;
epoch_start_end   = [0.25 1.049];% start and end of epoch, relative to trigger, in seconds
intertrial_trigger_num = 11;
blank_condition   = 10;

% Preprocess variables
var_threshold         = [.05 20]; % acceptable limits for variance in an epoch, relative to median of all epochs
bad_channel_threshold = 0.2;      % if more than 20% of epochs are bad for a channel, eliminate that channel
bad_epoch_threshold   = 0.2;      % if more than 20% of channels are bad for an epoch, eliminate that epoch
verbose               = false;

% condition names correspond to trigger numbers
condition_names  = {   ...
    'White Noise' ...
    'Binarized White Noise' ...
    'Pink Noise' ...
    'Brown Noise' ...
    'Gratings(0.36 cpd)' ...
    'Gratings(0.73 cpd)' ...
    'Gratings(1.45 cpd)' ...
    'Gratings(2.90 cpd)' ...
    'Plaid'...
    'Blank'};


% Where to find data?
project_pth    = '/Volumes/server/Projects/MEG/Gamma/Data';

% Type of data
data_pth       = '*_Gamma_*subj*';


% Find subject path
d = dir(fullfile(project_pth, data_pth));
%   restrict to directories
subj_pths = struct2cell(d);
isdir     = cell2mat(subj_pths(4,:));
subj_pths = subj_pths(1,isdir);

for subject = subjects

% Load denoised timeseries
data = load(fullfile(project_pth, subj_pths{subject}, 'processed',sprintf('s0%d_denoisedData.mat',subject+1)));
ts = data.denoisedts{1};
ts = permute(ts,[2,3,1]);

% Load bad channels and epoch matrices
badEpochs = data.bad_epochs;
badChannels = data.bad_channels;

% Get raw ts for triggers and then conditions again
raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject}, 'raw'), '*Gamma*');
trigger = meg_fix_triggers(raw_ts(:,trigger_channels));

% Get conditions
[~, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
% remove intertrial intervals
iti               = conditions == intertrial_trigger_num;
conditions        = conditions(~iti);
conditions_unique = unique(conditions);
num_conditions    = length(condition_names);

clear raw_ts

% Truncate design
conditions = conditions(~badEpochs);

% compute spectral data
t = (1:size(ts,1))/fs;
f = (0:length(t)-1)/max(t);

spectral_data = abs(fft(ts))/length(t)*2;

if nboot >= 1
    spectral_data_boots = zeros(size(ts,1), length(conditions_unique), size(ts,3), nboot);
else
    spectral_data_boots = zeros(size(ts,1), length(conditions_unique), size(ts,3));
end

% compute the mean amplitude spectrum for each electrode in each condition
fprintf('Computing bootstraps for each condition\n');
for ii = 1:length(conditions_unique)
    fprintf('Condition %d of %d\n', ii, length(conditions_unique)); drawnow;
    
    % Binary vector to identify epochs with this condition
    these_epochs = conditions == conditions_unique(ii);
    
    % spectral data, time points x epochs x channel
    these_data = spectral_data(:,these_epochs,:);
    
    if nboot > 1
        % reshape so that epochs are in rows (needed for bootstrp)
        these_data = permute(these_data, [2 1 3]);
        
        % log normalized mean of spectral power
        bootfun = @(x) squeeze(exp(nanmean(log(x),1)));
        
        % bootstat by definition is a matrix: nboot x (freq x channel)
        bootstat = bootstrp(nboot, bootfun, these_data);
        
        % reshape bootstat to 3D-array: nboot x freq x channel
        bootstat = reshape(bootstat, nboot, length(t), []);
        
        % spectral_data_boots is freq x condition x channel x boot
        spectral_data_boots(:,ii,:,:) = permute(bootstat,[2 3 1]);
        
    else
        spectral_data_boots(:,ii,:) = exp(nanmean(log(these_data),2));
    end
end
fprintf('Done!\n');

% Summarize bootstrapped spectral by mean and std over bootstraps
spectral_data_mean = mean(spectral_data_boots, 4);

%% Broadband and Gaussian Fit

% Convert the amplitude spectrum in each channel and each epoch into 2
% numbers, one for broadband and one for gamma

f_use4fit = f((f>=35 & f < 40) |(f > 40 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
f_sel=ismember(f,f_use4fit);
num_time_points = round((epoch_start_end(2)-epoch_start_end(1)+0.001)*fs);

num_channels = size(ts,3);
out_exp = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
gauss_f = NaN(num_channels,num_conditions, nboot);     % gaussian peak frequency
fit_f2  = NaN(num_conditions,num_time_points,num_channels, nboot); % fitted spectrum

warning off 'MATLAB:subsassigndimmismatch'

% For each channel, fit each condition separatley
fprintf('Fitting gamma and broadband values for each channel and each condition')
for cond = 1:num_conditions
    fprintf('Condition %d of %d\n', cond, length(conditions_unique)); drawnow;
    % Fit each channel separately
    for chan = 1:num_channels
        
        
        for bootnum = 1:nboot
            
            data_fit  = spectral_data_boots(:,cond,chan, bootnum);
            data_base = spectral_data_boots(:,blank_condition,chan, bootnum);
            
            % try/catch because bad channels / bad epochs were replaced by
            % NaNs, and NaNs will cause an error
            try
                [...
                    out_exp(chan, cond, bootnum), ...
                    w_pwr(chan, cond, bootnum), ...
                    w_gauss(chan, cond, bootnum),...
                    gauss_f(chan, cond, bootnum),...
                    fit_f2(cond,:, chan, bootnum)] = ...
                    gamma_fit_data(f,f_use4fit,data_base,data_fit);
            catch ME
                warning(ME.identifier, ME.message)
            end
        end
    end
end
fprintf('done!\n')

warning on 'MATLAB:subsassigndimmismatch'

% summarize bootstrapped fits
out_exp_mn = nanmean(out_exp,3);
w_pwr_mn   = nanmean(w_pwr,3);
w_gauss_mn = nanmean(w_gauss,3);
gauss_f_mn = nanmean(gauss_f,3);
fit_f2_mn  = nanmean(fit_f2,4);

out_exp_sd = nanstd(out_exp,[],3);
w_pwr_sd   = nanstd(w_pwr,[],3);
w_gauss_sd = nanstd(w_gauss,[],3);
gauss_f_sd = nanstd(gauss_f,[],3);
fit_f2_sd  = nanstd(fit_f2,[],4);

out_exp_md = nanmedian(out_exp,3);
w_pwr_md   = nanmedian(w_pwr,3);
w_gauss_md = nanmedian(w_gauss,3);
gauss_f_md = nanmedian(gauss_f,3);
fit_f2_md  = nanmedian(fit_f2,4);


fname = fullfile(project_pth, subj_pths{subject}, 'processed',sprintf('s0%d_denoisedData_bootstrapped100.mat',subject+1));
    parsave([fname '.mat'], 'out_exp', out_exp, 'w_pwr', w_pwr, ...
        'w_gauss', w_gauss, 'gauss_f', gauss_f,...
        'fit_f2', fit_f2, 'nboot', nboot);

end
%% Do some plotting

%% Power (mean condition)
fH = figure(996); clf, set(fH, 'name', 'Gaussian weight Before denoising')
for cond = 1:9
    subplot(3,3,cond)
    ft_plotOnMesh(to157chan(w_gauss_mn(:,cond)',~badChannels,0), condition_names{cond});
    set(gca, 'CLim', [0 .2])
end


fH = figure(997); clf, set(fH, 'name', 'Broadband weight Before denoising')
for cond = 1:9
    subplot(3,3,cond)
    ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)',~badChannels,0), condition_names{cond});
    set(gca, 'CLim', [0 2])
end

%% Weights (mean condition - baseline)
fH = figure(998); clf, set(fH, 'name', 'Gaussian weight')
for cond = 1:9
    subplot(3,3,cond)
    ft_plotOnMesh(to157chan(w_gauss_mn(:,cond)' - w_gauss_mn(:,num_conditions)', ~badChannels, 0), condition_names{cond});
    set(gca, 'CLim', [0 .2])
end


fH = figure(1001); clf
for cond = 1:9
    subplot(3,3,cond)
    ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
    set(gca, 'CLim', [-1 1] *.03)
end



