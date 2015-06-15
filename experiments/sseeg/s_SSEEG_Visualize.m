%% Visualize  SSEEG results after denoising

project_path          = '/Volumes/server/Projects/EEG/SSEEG/';
which_subject         = 'subj001';

%% Find the stored data
tmp = dir(fullfile(project_path, 'Data', 'Session*'));

subject_paths = struct2cell(tmp);
subject_paths = subject_paths(1,:);
which_path = ~cellfun(@isempty, strfind(subject_paths, which_subject));
subject_path = subject_paths{which_path};

% results of denoising broadband data
tmp = dir(fullfile(project_path, 'data', subject_path, 'processed', '*bb50.mat'));
if numel(tmp)>0,
    data_path_bb = fullfile(project_path, 'data', subject_path, 'processed', tmp(1).name);
end

% results of denoising stimulus locked data
tmp = dir(fullfile(project_path, 'data', subject_path, 'processed', '*sl40.mat'));
if numel(tmp)>0,
    data_path_sl = fullfile(project_path, 'data', subject_path, 'processed', tmp(1).name);
end

%% Load broadband results

bb=load(data_path_bb);

%% Plot broadband results

fH = figure(1); clf, set(fH, 'name', 'Denoised Broadband')
subplot(3,3,1)
data_to_plot = zeros(1, 128);
data_to_plot(~bb.badChannels) = bb.results.origmodel.r2;
plotOnEgi(data_to_plot), title('original R2'), colorbar; caxis([0 4]);
subplot(3,3,4)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2;
plotOnEgi(data_to_plot), title('final R2'), colorbar; caxis([0 4]);
subplot(3,3,7)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2 - bb.results.origmodel.r2;
plotOnEgi(data_to_plot, true), title('final R2 - original R2'), colorbar; caxis([0 4]);

a = [2 5 8];
cond = {'Full', 'Right', 'Left'};
for ii = 1:3
    subplot(3,3,a(ii))
    data_to_plot(~bb.badChannels) = bb.results.origmodel.beta_md(ii,:) ./ ...
        bb.results.origmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot); title(sprintf('SNR original %s ', cond{ii})); 
    colorbar; set(colorbar, 'Limits', [-2.5 2.5]);
    caxis([-2.5 2.5]);
    
    subplot(3,3,a(ii)+1);
    data_to_plot(~bb.badChannels) = bb.results.finalmodel.beta_md(ii,:) ./ ...
        bb.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii})), 
    colorbar; set(colorbar, 'Limits', [-2.5 2.5]);
    caxis([-2.5 2.5]);
end

%% Load stimulus locked results

sl = load(data_path_sl);

%% Plot stimulus locked results
fH = figure(2); clf, set(fH, 'Name', 'Stim-Locked Denoised');
    subplot(3,3,1)
data_to_plot = zeros(1, 128);
data_to_plot(~sl.badChannels) = sl.results.origmodel.r2;
plotOnEgi(data_to_plot), title('original R2'), colorbar; clim = get(subplot(3,3,1), 'CLim');

    subplot(3,3,4)
data_to_plot(~sl.badChannels) = sl.results.finalmodel.r2;
plotOnEgi(data_to_plot), title('final R2'), colorbar; set(subplot(3,3,4), 'CLim', clim);

    subplot(3,3,7)
data_to_plot(~sl.badChannels) = sl.results.finalmodel.r2 - sl.results.origmodel.r2;
plotOnEgi(data_to_plot), title('final R2 - original R2'), colorbar; set(subplot(3,3,7), 'CLim', clim);

a = [2 5 8];
cond = {'Full', 'Right', 'Left'};
for ii = 1:3
    subplot(3,3,a(ii))
    data_to_plot(~sl.badChannels) = sl.results.origmodel.beta_md(ii,:) ./ ...
        sl.results.origmodel.beta_se(ii,:);        
    plotOnEgi(data_to_plot), title(sprintf('SNR original %s ', cond{ii})); 
    colorbar; set(colorbar, 'Limits', [0 10]);
    caxis([0 10]);
    
    subplot(3,3,a(ii)+1);
    data_to_plot(~sl.badChannels) = sl.results.finalmodel.beta_md(ii,:) ./ ...
        sl.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii}));
    colorbar; set(colorbar, 'Limits', [0 10]);
    caxis([0 10]);
end

%% Visualize the noise pool

noise_pool = zeros(1,128);
noise_pool(results_60_noise.noisepool) = true;
figure; plotOnEgi(noise_pool); title('Noise pool');

%% Headplot of difference in SNR between left and right conditions

cond_diff_headplot(bb, 'left_right', 'finalmodel', 'SNR');

%% Plot spectra of non-denoised versus denoised timeseries

ts_orig = sensorData;
ts_den = permute(sl.denoisedts{1}, [2 3 1]);
t = size(ts,1);
num_epoch_time_pts = t;
data_channels = 1:128;
channels_to_plot = [84 128] ;
produce_figures = 1;

addpath(genpath('/Volumes/server/Projects/MEG/SSMEG/code/'));
freq                    = (0:num_epoch_time_pts-1)/(num_epoch_time_pts/1000); 
on_full_orig            = ts_orig(:, find(conditions == 1), :); 
on_full_den             = ts_den(:, find(conditions == 1), :); 
ft_on_epoched_full_orig = fft(on_full_orig) / length(t)*2;   
ft_on_epoched_full_den  = fft(on_full_den) / length(t)*2;   
amps_on_full_orig       = abs(ft_on_epoched_full_orig);    clear ft_on_epoched_full;
amps_on_full_den        = abs(ft_on_epoched_full_den);    clear ft_on_epoched_full;

    figure(206);
    clf; set(gcf, 'Color', 'w')
    set(gca, 'FontSize', 20, 'ColorOrder', jet(length(channels_to_plot)))
    hold all;
    pre_84 = plot(freq, squeeze(nanmedian(amps_on_full_orig(:,:,channels_to_plot(1)), 2)), 'LineWidth', 2, 'Color', 'r');
    post_84 = plot(freq, squeeze(nanmedian(amps_on_full_den(:,:,channels_to_plot(1)), 2)), 'LineWidth', 2, 'Color', 'g');
    pre_128 = plot(freq, squeeze(nanmedian(amps_on_full_orig(:,:,channels_to_plot(2)), 2)), 'LineWidth', 2, 'Color', 'b');
    post_128 = plot(freq, squeeze(nanmedian(amps_on_full_den(:,:,channels_to_plot(2)), 2)), 'LineWidth', 2, 'Color', 'c');
    xlim([0 85])
    yl = get(gca, 'YLim');
    for ii = 1:7; plot(ii*freq(1)*[12 12], yl, 'k-'); end
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (microvolts)')
    legend('channel 84 pre', 'channel 84 post','channel 128 pre', 'channel 128 post');
    title('Subj001 channel 84 and 128 spectra before and after stim-locked denoising')
    get_MEG_axes('True'); 
    sub_head = axes('position', [.15 .2 .2 .2]);
    temp = zeros(1,128); temp(channels_to_plot) = 1;
    plotOnEgi(temp);

%% Visually compare EEG data from before and after denoising

visual_channels = 55:95;
ts_cat = [];
denoised_ts_cat = [];

for ii = 1:72
    % ts_cat          = cat(2, ts_cat, sensorData(visual_channels, :, ii));
    denoised_ts_cat = cat(2, denoised_ts_cat, sl.denoisedts{1}(visual_channels, :, ii));
end

% figure(13); plot(ts_cat'); title('Original'); hold all;
figure(10); plot(denoised_ts_cat'); title('Denoised electrode data');
xlabel('frames (ms)')
ylabel('Voltage (microvolts?)')
