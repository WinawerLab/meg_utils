% s_GAMMA_MEG_analysis

% Analyze and visualize data from MEG Gamma experiments. Subjects saw
% several kinds of stimuli, including gratings of various spatial
% frequencies, plaids, and noise patterns (phase-randomized with 1/f^n
% spectral power distrubutions)
%
% Stimuli were static, and on the screen either for 500 ms (subjects 1-3)
% or 1000 ms (subjects 4-6) with 500 ms ISI.
%
% Spectral data from each channel for each stimulus type are modeled as a
% mixture of a line and gaussian in log power / log frequency (meaning, a
% power law and a narrowband response)
%
% TODO:
%   1. Summarize the responses from each subject, for example as images on
%   meshes showing the gaussian response and the broadband response for
%   each separate type of image, as well as for contrasts of several
%   stimulus categories (e.g., grating v noise, or grating v blank)
%
%   2. Denoise with environmental noise channels
%
%   3. MAYBE: summarize the EVOKED response from each stimulus
%
%   4. Make a POSTER!
% 
%   5. MAYBE: a grand average where you average the mesh images across
%   subjects
%
%   6. MAYBE: source localize the signals

% Analysis options
%% Set analysis variables
project_pth                     = '/Volumes/server/Projects/MEG/Gamma/Data';

% data to be analysed
data_pth                      = '*_Gamma_*subj*';

data_channels                 = 1:157;
environmental_channels        = 158:160;
trigger_channels              = 161:164;

denoise_with_nonphys_channels = true;        % Regress out time series from 3 nuissance channels
remove_bad_epochs             = true;        % Remove epochs whose variance exceeds some threshold
remove_bad_channels           = true;        % Remove channels whose median sd is outside some range

produce_figures               = true;        % If you want figures in case of debugging, set to true

denoise_via_pca               = false;       % Do you want to use megdenoise?

fs                            = 1000;        % sample rate
epoch_start_end               = [0.550 1.05];% start and end of epoch, relative to trigger, in seconds

intertrial_trigger_num        = 11;          % the MEG trigger value that corresponds to the intertrial interval

save_images                   = false;

% condition names correspond to trigger numbers
condition_names               = {   ...
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

which_data_sets_to_analyze = 1;
blank_condition = strcmpi(condition_names, 'blank');
%% Add paths

%change server-1 back to server
meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip', 'yokogawa_defaults')

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);
%% Loops over datasets
for subject_num = which_data_sets_to_analyze
    
    save_pth = fullfile(project_pth, 'Images', subj_pths{subject_num});
    if ~exist(save_pth, 'dir'), mkdir(save_pth); end
    
    % --------------------------------------------------------------------
    % ------------------ PREPROCESS THE DATA -----------------------------
    % --------------------------------------------------------------------
    %% Load data (SLOW)
    raw_ts = meg_load_sqd_data(fullfile(project_pth, subj_pths{subject_num}, 'raw'), '*Gamma*');
    
    %% Extract triggers
    trigger = meg_fix_triggers(raw_ts(:,trigger_channels));
    
    %% Make epochs
    [ts, conditions]  = meg_make_epochs(raw_ts, trigger, epoch_start_end, fs);
    % remove intertrial intervals
    iti               = conditions == intertrial_trigger_num;
    ts                = ts(:,~iti, :);
    conditions        = conditions(~iti);
    conditions_unique = unique(conditions);
    num_conditions    = length(condition_names);
    
    %% Find bad epochs
    if remove_bad_epochs
        
        % This identifies any epochs whos variance is outside some multiple of the
        % grand variance
        bad_epochs = meg_find_bad_epochs(ts(:,:,data_channels), [.05 20]);
        
        % any epoch in which more than 10% of channels were bad should be removed
        % entirely
        epochs_to_remove = mean(bad_epochs,2)>.1;
        
        % once we remove 'epochs_to_remove', check whether any channels have more
        % than 10% bad epochs, and we will remove these
        channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.1;
        
        bad_epochs(epochs_to_remove,:) = 1;
        bad_epochs(:,channels_to_remove) = 1;
        
        figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')
        
        ts = meg_remove_bad_epochs(bad_epochs, ts);
    end
    
    
    %% Denoise data by regressing out nuissance channel time series
    
    % TODO: check whether this runs with NaNs in ts
    
    % Denoise data with 3 noise channels
    if denoise_with_nonphys_channels
        if exist('./denoised_with_nuissance_data.mat', 'file')
            load(fullfile(data_pth{subject_num},'denoised_with_nuissance_data.mat'));
        else fprintf('Loading data.. This may take a couple of seconds\n');
            ts = meg_environmental_denoising(ts, environmental_channels,...
                data_channels, produce_figures);
        end
    end
    
    
    % --------------------------------------------------------------------
    % ------------------ ANALYZE THE PREPROCESSED DATA -------------------
    % --------------------------------------------------------------------
    %% Spectral analysis
    
    % compute spectral data
    t = (1:size(ts,1))/fs;
    f = (0:length(t)-1)/max(t);
    nboot = 3; % number of bootstrap samples
    spectral_data = abs(fft(ts))/length(t)*2;
    spectral_data_boots = zeros(size(ts,1), length(conditions_unique), length(data_channels), nboot);
    
    % compute the mean amplitude spectrum for each electrode in each condition
    fprintf('Computing bootstraps for each condition');
    for ii = 1:length(conditions_unique)
        fprintf('.'); drawnow;
        
        % Binary vector to identify epochs with this condition
        these_epochs = conditions == conditions_unique(ii);
        
        % spectral data, time points x epochs x channel
        these_data = spectral_data(:,these_epochs,data_channels);
        
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
        
    end
    fprintf('Done!\n');
    
    % Summarize bootstrapped spectral by mean and std over bootstraps
    spectral_data_mean = mean(spectral_data_boots, 4);
    spectral_data_std  =  std(spectral_data_boots, [], 4);
    spectral_data_snr   = spectral_data_mean./spectral_data_std;
    
    %% Broadband and Gaussian Fit
    
    % Convert the amplitude spectrum in each channel and each epoch into 2
    % numbers, one for broadband and one for gamma
    
    f_use4fit = f((f>=35 & f <= 57) | (f>=65 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
    f_sel=ismember(f,f_use4fit);
    
    num_channels = length(data_channels);
    out_exp = NaN(num_channels,num_conditions, nboot);     % slope of spectrum in log/log space
    w_pwr   = NaN(num_channels,num_conditions, nboot);     % broadband power
    w_gauss = NaN(num_channels,num_conditions, nboot);     % gaussian height
    gauss_f = NaN(num_channels,num_conditions, nboot);     % gaussian peak frequency
    fit_f2  = NaN(num_conditions,1000,num_channels, nboot); % fitted spectrum
    
    % For each channel, fit each condition separatley
    fprintf('Fitting gamma and broadband values for each channel and each condition')
    for cond = 1:num_conditions
        fprintf('.'); drawnow;
        % Fit each channel separately
        for chan = data_channels
            
            
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
    
    % summarize bootstrapped fits
    out_exp_mn = mean(out_exp,3);
    w_pwr_mn   = mean(w_pwr,3);
    w_gauss_mn = mean(w_gauss,3);
    gauss_f_mn = mean(gauss_f,3);
    fit_f2_mn  = mean(fit_f2,4);
    
    out_exp_sd = std(out_exp,[],3);
    w_pwr_sd   = std(w_pwr,[],3);
    w_gauss_sd = std(w_gauss,[],3);
    gauss_f_sd = std(gauss_f,[],3);
    fit_f2_sd  = std(fit_f2,[],4);

    
    %% Plot Gaussian fits
    line_width = 2; % line width for
    for chan = data_channels
        fH = figure(10); clf, set(gcf, 'Position', [100 100 800 800], 'Color', 'w')
        
        data_base = spectral_data_mean(:,num_conditions,chan);
        for cond = 1:num_conditions-1
            data_fit = spectral_data_mean(:,cond,chan,:);
            subplot(3,3,cond)
            % plot fit
            plot(f,10.^(fit_f2_mn(cond,:, chan)),'Color','g','LineWidth',line_width)
            hold on;
            
            % plot baseline data
            plot(f(f_sel),data_base(f_sel),'k--','LineWidth',line_width)
            % plot stimulus data
            plot(f(f_sel),data_fit(f_sel),'-','LineWidth',line_width)
            
            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                'XLim', [min(f_use4fit) max(f_use4fit)], ...
                'YLim', 10.^[0.3 1.5] )
            title(sprintf('Subject %d, Channel %d, %s', subject_num, chan, condition_names{cond}))            
        end
        if save_images,
            
            hgexport(fH, fullfile(save_pth, sprintf('Spectra_Chan%03d.eps', chan)));
        else
            pause(0.1)
        end
    end
    
    %% Plot spectra
    fH = figure(2); clf,  set(gcf, 'Position', [100 100 1000 400], 'Color', 'w')
    
    colors = zeros(num_conditions,3);
    colors(1:4,:) = ([1 1 1]'*[.2 .4 .6 .8])';
    colors(5:9,:) = hsv(5);
    colors(10,:)  = [0 0 0];
    yl            =  10.^([0.5 1.6]);
    xl            = [30 200];
    
    for chan = data_channels
        
        % Plot Data -----------------------------------------------------------
        subplot(1,3,1); cla; hold on
        for cond = 1:num_conditions
            plot(f(f_sel), squeeze(spectral_data_mean(f_sel,cond,chan)), 'Color', colors(cond,:));
        end
        set(gca, 'YScale', 'log','XScale', 'log', 'YLim', yl, 'XLim', xl,  ...
            'Color', [1 1 1], 'XGrid', 'on', ...
            'XTick', [10 60 100 200]);
        title(sprintf('Data from channel %d, subject %d', chan, subject_num))
        
        % Plot Fits -----------------------------------------------------------
        subplot(1,3,2); cla; hold on
        for cond = 1:num_conditions
            plot(f, squeeze(10.^(fit_f2_mn(cond,:, chan))), 'Color', colors(cond,:));
        end
        set(gca, 'YScale', 'log','XScale', 'log', 'YLim', yl, 'XLim', xl, ...
            'Color', [1 1 1], 'XGrid', 'on', ...
            'XTick', [10 60 100 200]);
        title(sprintf('Fits to channel %d, subject %d', chan, subject_num))
        
        % Legend  -----------------------------------------------------------
        subplot(1,3,3); cla
        
        set(gca, 'ColorOrder', colors); hold all
        plot(zeros(10,10), zeros(10,10), '-')
        box off; axis off;
        legend(condition_names)
        drawnow;
        if save_images
            hgexport(fH, fullfile(save_pth, sprintf('Spectral_fits_Chan%03d.eps', chan)));
        end
        
    end
    %axis tight
    
    %% Mesh visualization of model fits
    
    % TODO: threshold maps by significance: w_gauss_mn./w_gauss_sd>2
    
    
    fH = figure(998); clf, set(fH, 'name', 'Gaussian weight')
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(w_gauss_mn(:,cond)', condition_names{cond});
        set(gca, 'CLim', [0 .2])
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Gaussian.eps'));
    end
    
    fH = figure(999); clf
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)', condition_names{cond});
        set(gca, 'CLim', [-1 1] *.03)
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Broadband.eps'));
    end
    
    
    fH = figure(1000); clf
    subplot(2,2,1)
    ft_plotOnMesh((w_gauss * [0 0 0 0 1 1 1 1 0 -4]')', 'Gamma power, All gratings minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,2)
    ft_plotOnMesh((w_gauss * [1 1 1 1 0 0 0 0 0 -4]')', 'Gamma power, All noise minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,3)
    ft_plotOnMesh((w_gauss * [-1 -1 -1 -1 1 1 1 1 0 0]')', 'Gamma power, All gratings minus all noise');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,4)
    ft_plotOnMesh((w_pwr * [1 1 1 1 1 1 1 1 1 -9]')', 'Broadband, All stimuli minus baseline');
    set(gca, 'CLim', .2 * [-1 1])
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Gamma_Gratings_M_Baseline.eps'));
    end
    
    
end



