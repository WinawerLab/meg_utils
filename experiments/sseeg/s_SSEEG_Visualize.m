%% Visualize  SSEEG results after denoising

project_path          = '/Volumes/server/Projects/EEG/SSEEG/';
which_subject         = 'subj004';

%% Find the stored data
tmp = dir(fullfile(project_path, 'Data', 'Session*'));

subject_paths = struct2cell(tmp);
subject_paths = subject_paths(1,:);
which_path = ~cellfun(@isempty, strfind(subject_paths, which_subject));
subject_path = subject_paths{which_path};

% results of denoising broadband data
tmp = dir(fullfile(project_path, 'data', subject_path, 'processed', '*bb.mat'));
if numel(tmp)>0,
    data_path_bb = fullfile(project_path, 'data', subject_path, 'processed', tmp(1).name);
end

% results of denoising stimulus locked data
tmp = dir(fullfile(project_path, 'data', subject_path, 'processed', '*sl.mat'));
if numel(tmp)>0,
    data_path_sl = fullfile(project_path, 'data', subject_path, 'processed', tmp(1).name);
end

%% Plot broadband results

bb=load(data_path_bb);

fH = figure(1); clf, set(fH, 'name', 'Denoised Broadband')
subplot(3,3,1)
data_to_plot = zeros(1, 128);
data_to_plot(~bb.badChannels) = bb.results.origmodel.r2;
plotOnEgi(data_to_plot), title('original R2'), colorbar; clim = get(subplot(3,3,1), 'CLim');
subplot(3,3,4)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2;
plotOnEgi(data_to_plot), title('final R2'), colorbar; set(subplot(3,3,4), 'CLim', clim);
subplot(3,3,7)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2 - bb.results.origmodel.r2;
plotOnEgi(data_to_plot), title('final R2 - original R2'), colorbar; set(subplot(3,3,7), 'CLim', clim);

a = [2 5 8];
cond = {'Full', 'Right', 'Left'};
for ii = 1:3
    subplot(3,3,a(ii))
    data_to_plot(~bb.badChannels) = bb.results.origmodel.beta_md(ii,:) ./ ...
        bb.results.origmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR original %s ', cond{ii})), colorbar;
    clim = get(subplot(3,3,a(1)), 'CLim'); clim(2) = -clim(1);
    set(subplot(3,3,a(ii)), 'CLim', clim);
    
    subplot(3,3,a(ii)+1);
    data_to_plot(~bb.badChannels) = bb.results.finalmodel.beta_md(ii,:) ./ ...
        bb.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii})), colorbar;
    set(subplot(3,3,a(ii)+1), 'CLim', clim);
end

%% Plot stimulus locked results

sl = load(data_path_sl);

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
    colorbar; clim = get(subplot(3,3,a(1)), 'CLim'); clim(2) = -clim(1);
    set(subplot(3,3,a(ii)), 'CLim', clim);
    
    subplot(3,3,a(ii)+1);
    data_to_plot(~sl.badChannels) = sl.results.finalmodel.beta_md(ii,:) ./ ...
        sl.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii}));
    colorbar; set(subplot(3,3,a(ii)+1), 'CLim', clim);
end

%% Visualize the noise pool

noise_pool = zeros(1,128);
noise_pool(bb.results.noisepool) = true;
figure; plotOnEgi(noise_pool); title('Noise pool');

%% Visualize
% sseegMakePrePostHeadplot(project_path,session_name,session_prefix,true)

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


%% Plot spectra of non-denoised timeseries

denoisedts = reshape(denoisedts{1},2,3,1);
t = size(denoisedts,1);
num_epoch_time_pts = t;
data_channels = 1:128;
channels_to_plot = 70:75;
produce_figures = 1;


[amps_on_full,amps_on_right,amps_on_left, amps_off_full,amps_off_right, ...
    amps_off_left] = sseeg_fourier(t, num_epoch_time_pts, denoisedts, conditions, ...
       data_channels, channels_to_plot, produce_figures);

