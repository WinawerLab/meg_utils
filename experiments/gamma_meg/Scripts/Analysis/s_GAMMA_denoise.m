%% denoiseGamma scripts
%
% DESCRIPTION: Wrapper script to denoise  MEG visual gamma data sets.


% ------------------------------------------------------------------------
% ----------------- Check options & define variables----------------------
% ------------------------------------------------------------------------

% Where to find data?
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
ft_pth                        = '/Volumes/server/Projects/MEG/code/fieldtrip';

addpath(genpath('~/matlab/git/denoiseproject'))

% Type of data
data_pth                      = '*_Gamma_*subj*';

% Subject number to analyze
subjects                      = 17; % 99 means synthetic data (remember to set denoise_with_nonphys_channels to false)

% preprocessing parameters (see dfdPreprocessData)
var_threshold                 = [0.05 20];
bad_channel_threshold         = 0.2;
bad_epoch_threshold           = 0.2;
data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;
epoch_start_end               = [0.050 1.049];% start and end of epoch, relative to trigger, in seconds
fs                            = 1000;        % sample rate
intertrial_trigger_num        = 14;          % the MEG trigger value that corresponds to the intertrial interval
blank_condition               = 13;          % the MEG trigger value that corresponds to trials with zero contrast
verbose                       = true;
denoise_with_nonphys_channels = true;
save_data                     = true;


% Find subject path
d = dir(fullfile(project_pth, data_pth));
%   restrict to directories
subj_pths = struct2cell(d);
isdir     = cell2mat(subj_pths(4,:));
subj_pths = subj_pths(1,isdir);

% ------------------------------------------------------------------------
% --------------------- Load & preprocess data ---------------------------
% ------------------------------------------------------------------------
for subject = subjects    

    %% check paths
    if isempty(which('sqdread')),
        meg_add_fieldtrip_paths(ft_pth,{'yokogawa', 'sqdproject'})
    end
    
    
    if subject == 99
        [ts, conditions] = gamma_synthetize_validation_data();
    else
        %% Load data (SLOW)
        raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject}, 'raw'), '*Gamma*');
        
        %% Extract triggers
        trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
        
        %% Make epochs
        %   ts is [time points by epochs x channels]
        [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
        % remove intertrial intervals
        iti               = conditions == intertrial_trigger_num;
        ts                = ts(:,~iti, :);
        conditions        = conditions(~iti);
    end
    
    condition_names = gamma_get_condition_names(subject);    
    conditions(conditions==blank_condition) = 0;
    
    %% Remove bad epochs and bad channels
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels
        ts = meg_environmental_denoising(ts, environmental_channels,data_channels);
    end
    
    [ts(:,:,data_channels), bad_channels, bad_epochs] = meg_preprocess_data(ts(:,:,data_channels), ...
        var_threshold, bad_channel_threshold, bad_epoch_threshold, 'meg160xyz', verbose);
    
    % clip out the bad channels and bad epochs from the data, and remove the
    % bad epochs from the conditions vector
    ts = ts(:,~bad_epochs,~bad_channels);
    conditions = conditions(~bad_epochs);
    
    % convert the conditions vector to a binary design matrix
    design_mtrx = conditions2design(conditions);
    
    %% -----------------------------------------------------------------------
    % -------------------------- Denoise data --------------------------------
    % ------------------------------------------------------------------------
    
    % Define denoise Parameters (see denoisedata.m)
    opt.pcchoose        = 1.05;
    opt.npcs2try        = 10;
    opt.resampling      = {'xval','xval'}; % could be {'boot' 'boot'};
    opt.pcselmethod     = 'r2';            % could be 'snr';
    if verbose; opt.verbose = true; end
    
    % time and frequencies for full-length epochs (prior to truncating)
    evoked_cutoff       = .250; % use only time points >= evoked_cutoff
    t                   = epoch_start_end(1):1/fs:epoch_start_end(2);
    t_clipped           = t(t >= evoked_cutoff);
    f                   = (0:length(t_clipped)-1)*fs/length(t_clipped);
%     keep_frequencies    = @(x) x((x>=35 & x < 40) |(x > 40 & x <= 57) | ...
%         (x>=65 & x <= 115) | (x>=126 & x<= 175) | (x>=186 & x<= 200));

     keep_frequencies = @(x) x((x>=35 & x <= 57) | (x>=63 & x <= 115) | (x>=126 & x <= 175) | (x>=186 & x <= 200));
    
    % preprocess data by clipping 1st n time points (to eliminate evoked
    % response) and by removing all frequencies that will not be used to fit
    % gamma and broadband later on.
    opt.preprocessfun   = @(x)gammapreprocess(x, t, f, evoked_cutoff, keep_frequencies(f));
    
    % Evoked function is used to identify the noise pool. It finds the peak
    % evoked signal in each epoch for each channel.
    evokedfun           = @(x)getevoked(x, fs, design_mtrx, [-30 30]); % function handle to determine noise pool
    
    % Eval function that extracts broadband
    evalfun             = @(x)getbroadband(x,keep_frequencies, fs);  % function handle to compute broadband
    
    % Permute sensorData for denoising
    %   [time points by epochs x channels] => [channels x time points x epochs]
    ts = permute(ts, [3 1 2]);
    
    % Denoise for broadband analysis
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design_mtrx,ts,evokedfun,evalfun,opt);
    
    if save_data
        save(fullfile(project_pth, sprintf('s0%d_denoisedData.mat',subject+1)),'results','evalout','bad_channels','bad_epochs','denoisedts','opt')
    end
    
    %% Look at results 
    %
    % Note that the broadband plotted in these meshes are probably a
    % combination of Gamma oscillatory responses and Broadband responses.
    % For a good understanding of the data, run
    % s_Gamma_fit_data_afte_denoising.m. These figures plotted below are
    % just a sanity check.
    figure, ft_plotOnMesh(to157chan(results.noisepool, ~bad_channels, 0), ...
        'Noise Pool', [],'2d', 'interpolation', 'nearest');
    
%     figure, ft_plotOnMesh(to157chan(results.finalmodel.r2, ~bad_channels, 0), ...
%         'Denoised Broadband R2',  [], [],  'CLim', [0 5]);
%     
%     figure, ft_plotOnMesh(to157chan(results.origmodel.r2, ~bad_channels, 0), ...
%         'Original Broadband R2', [], [],  'CLim', [0 5]);
%     
%     figure, ft_plotOnMesh(to157chan(evalout(1).r2, ~bad_channels, 0), ...
%         'Broadband R2 0 PCs', [], [], 'CLim', [0 5]);
%     
%     
%     snr.final = results.finalmodel.beta_md  ./ results.finalmodel.beta_se;
%     figure, set(gcf, 'Name', 'Denoised Broadband SNR');
%     for ii = 1:9
%         subplot(3,3,ii)
%         ft_plotOnMesh(to157chan(snr.final(ii,:), ~bad_channels, 0), ...
%             condition_names{ii},  [], [], 'CLim', [-3 3]);
%     end
%     
%     
%     snr.orig = results.origmodel.beta_md  ./ results.origmodel.beta_se;
%     figure, set(gcf, 'Name', 'Original Broadband SNR');
%     for ii = 1:9
%         subplot(3,3,ii)
%         ft_plotOnMesh(to157chan(snr.orig(ii,:), ~bad_channels, 0), ...
%             condition_names{ii},  [], [], 'CLim', [-3 3]);
%     end
%     
%     figure, set(gcf, 'Name', 'Broadband SNR (final - orig)');
%     for ii = 1:9
%         subplot(3,3,ii)
%         ft_plotOnMesh(to157chan(snr.final(ii,:) - snr.orig(ii,:), ~bad_channels, 0), ...
%             condition_names{ii},  [], [], 'CLim', [-3 3]);
%     end
    
end
