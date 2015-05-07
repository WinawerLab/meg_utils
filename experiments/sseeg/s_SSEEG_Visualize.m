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


fH = figure(10); clf, set(fH, 'name', 'Denoised Broadband 3pcs')
subplot(3,3,1)
data_to_plot = zeros(1, 128);
data_to_plot(~bb.badChannels) = bb.results.origmodel.r2;
plotOnEgi(data_to_plot), title('original R2'), colorbar
subplot(3,3,4)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2;
plotOnEgi(data_to_plot), title('final R2'), colorbar
subplot(3,3,7)
data_to_plot(~bb.badChannels) = bb.results.finalmodel.r2 - bb.results.origmodel.r2;
plotOnEgi(data_to_plot), title('final R2 - original R2'), colorbar

a = [2 5 8];
cond = {'Full', 'Right', 'Left'};
for ii = 1:3
    subplot(3,3,a(ii))
    data_to_plot(~bb.badChannels) = bb.results.origmodel.beta_md(ii,:) ./ ...
        bb.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR original %s ', cond{ii})), colorbar;
    subplot(3,3,a(ii)+1);
    data_to_plot(~bb.badChannels) = bb.results.finalmodel.beta_md(ii,:) ./ ...
        bb.results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii})), colorbar;
end

%% Plot stimulus locked results

fH = figure(3); clf, set(fH, 'Name', 'Stim-Locked Denoised 6pcs')
subplot(3,3,1)
data_to_plot = zeros(1, 128);
data_to_plot(~badChannels) = results.origmodel.r2;
plotOnEgi(data_to_plot), title('original R2'), colorbar; clim = get(subplot(3,3,1), 'CLim');
subplot(3,3,4)
data_to_plot(~badChannels) = results.finalmodel.r2;
plotOnEgi(data_to_plot), title('final R2'), colorbar; set(subplot(3,3,4), 'CLim', clim);
subplot(3,3,7)
data_to_plot(~badChannels) = results.finalmodel.r2 - results.origmodel.r2;
plotOnEgi(data_to_plot), title('final R2 - original R2'), colorbar; set(subplot(3,3,4), 'CLim', clim);

a = [2 5 8];
cond = {'Full', 'Right', 'Left'};
for ii = 1:3
    subplot(3,3,a(ii))
    data_to_plot(~badChannels) = results.origmodel.beta_md(ii,:) ./ ...
        results.origmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR original %s ', cond{ii})), colorbar; clim = get(subplot(3,3,a(ii)), 'CLim');
    subplot(3,3,a(ii)+1);
    data_to_plot(~badChannels) = results.finalmodel.beta_md(ii,:) ./ ...
        results.finalmodel.beta_se(ii,:);
    plotOnEgi(data_to_plot), title(sprintf('SNR final %s ', cond{ii})), colorbar; set(subplot(3,3,a(ii)+1), 'CLim', clim);
end

%% Visualize the noise pool

noise_pool = zeros(1,128);
noise_pool(results.noisepool) = true;
figure; plotOnEgi(noise_pool); title('Noise pool');

%% Visualize
% sseegMakePrePostHeadplot(project_path,session_name,session_prefix,true)

%% Visually compare EEG data from before and after denoising

visual_channels = 55:95;
ts_cat = [];
denoised_ts_cat = [];

for ii = 1:24
    ts_cat          = cat(2, ts_cat, sensorData(visual_channels, :, ii));
    denoised_ts_cat = cat(2, denoised_ts_cat, denoisedts{1}(visual_channels, :, ii));
end

figure(13); plot(ts_cat'); title('Original');
figure(14); plot(denoised_ts_cat'); title('Denoised');


