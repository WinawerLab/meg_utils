%% GAMMA MEG make across subjects headplot

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';

fs                            = 1000;
intertrial_trigger_num        = 10;
session_num                   = [9:11,13];
save_images                   = false;
save_pth                      = fullfile(project_pth, 'Images');
use_denoised_data             = true;
if use_denoised_data;
     delete_me = 3; 
else delete_me = 0; 
end;

addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'))

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

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

AllResults   = [];
bad_channels = [];
snr_w_pwr    = [];
snr_w_gauss  = [];

for subject_num = session_num
    

    
    %% Load Data
    load_pth                = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    if use_denoised_data; d =  dir(fullfile(load_pth, '*denoisedData_bootstrapped100_2*'));  % Check whether d(1) or d(2) is used, for positive or both sides modelfit
    else                  d =  dir(fullfile(load_pth, '*_bootstrappedData.mat')); end
        
    AllResults{subject_num-delete_me} = load(fullfile(load_pth, d(1).name));
    w_gauss                 = AllResults{subject_num -delete_me}.w_gauss;
    w_pwr                   = AllResults{subject_num -delete_me}.w_pwr;
    
    % Load denoisedData to get bad channels.
    if use_denoised_data
        d                       =  dir(fullfile(load_pth, '*denoisedData.mat'));
        denoisedData            = load(fullfile(load_pth, d(1).name));
        bad_channels{subject_num-delete_me} = denoisedData.bad_channels;
    else bad_channels{subject_num} = zeros(1,157);
    end
    
    
    
    num_conditions          = size(AllResults{subject_num-delete_me}.out_exp,2);
    num_channels            = size(AllResults{subject_num-delete_me}.out_exp,1);
    
    %% Calculating SNR contrasts
    
    summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
    
    contrasts = [...
        1 0 0 0 0 0 0 0 0 -1; ...    % white noise - baseline
        0 1 0 0 0 0 0 0 0 -1; ...    % binarized white noise - baseline
        0 0 1 0 0 0 0 0 0 -1; ...    % pink noise - baseline
        0 0 0 1 0 0 0 0 0 -1; ...    % brown noise - baseline
        0 0 0 0 1 0 0 0 0 -1; ...    % 0.36cpd gratings - baseline
        0 0 0 0 0 1 0 0 0 -1; ...    % 0.73cpd gratings - baseline
        0 0 0 0 0 0 1 0 0 -1; ...    % 1.46cpd gratings - baseline
        0 0 0 0 0 0 0 1 0 -1; ...    % 2.90cpd gratings - baseline
        0 0 0 0 0 0 0 0 1 -1; ...    % plaid - baseline
        1 1 1 1 0 0 0 0 0 -4; ...    % noise - baseline
        0 0 0 0 1 1 1 1 0 -4; ...    % gratings - baseline
        1 1 1 1 -1 -1 -1 -1 0 0; ... % noise - gratings
        -1 -1 -1 -1 1 1 1 1 0 0; ... % gratings - noise
        ];
    contrastnames = {
        'white noise - baseline'...
        'binwn - baseline'...
        'pink noise - baseline'...
        'brown noise - baseline'...
        '0.36cpd gratings - baseline'...
        '0.73cpd gratings - baseline'...
        '1.46cpd gratings - baseline'...
        '2.90cpd gratings - baseline'...
        'plaid - baseline' ...
        'noise - baseline'...
        'gratings - baseline'...
        'noise - gratings'...
        'gratings - noise'...
        };
    
    
    % ensure each condition is weighted proportionatly in each contrast
    contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
    
    num_contrasts = size(contrasts,1);
    
    % compute SNR
    

        
    tmp_data = permute(w_pwr, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, AllResults{subject_num-delete_me}.nboot);
    snr_w_pwr{subject_num-delete_me} = summary_stat(tmp)';
    
    tmp_data = permute(w_gauss, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, AllResults{subject_num-delete_me}.nboot);
    snr_w_gauss{subject_num-delete_me}  = summary_stat(tmp)';
    
end

snr_w_gauss_157 = [];
snr_w_pwr_157   = [];
for ii = 1:num_contrasts
    for jj = session_num-delete_me
        snr_w_gauss_157{ii,jj} = to157chan(snr_w_gauss{jj}(:,ii)',~bad_channels{jj}, 'nans');
        snr_w_pwr_157{ii,jj} = to157chan(snr_w_pwr{jj}(:,ii)',~bad_channels{jj}, 'nans');      
    end
end

snr_w_gauss_157 = reshape(catcell(1,snr_w_gauss_157),[num_contrasts,numel(session_num),157]);
snr_w_pwr_157   = reshape(catcell(1,snr_w_pwr_157),[num_contrasts,numel(session_num),157]);


%% SNR Mesh per stimulus type Gamma
threshold = 0;%3;
% gaussian weight for each stimuli
fH = figure(1005); clf, set(fH, 'name', 'Gaussian weight' )

for c = 1:9
    subplot(3,3,c)
    data_to_plot = squeeze(nanmean(snr_w_gauss_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]*2)
end

data = {};
data.snr_w_gauss   = snr_w_gauss;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};

hgexport(fH, fullfile(save_pth,'Mesh_stimulus_type_gamma_SNR_across_subjects_BinarizedDenoised_4'));
set(gcf, 'UserData', data);
saveas(fH, fullfile(save_pth,'Mesh_stimulus_type_gamma_SNR_across_subjects_BinarizedDenoised_4'), 'fig');

%% SNR Mesh difference Noise and Gratings
threshold = 0;%3;
% gaussian weight for each stimuli
fH = figure(455); clf, set(fH, 'name', 'Gaussian weight')
%for c = 1:12
for c = 10:13
    subplot(2,2,c-9)
    data_to_plot = squeeze(nanmean(snr_w_gauss_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]* 2)
end

data = {};
data.snr_w_gauss   = snr_w_gauss;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};

hgexport(fH, fullfile(save_pth,'Mesh_difGratNoise_gamma_SNR_across_subjects_BinarizedDenoised_4'));
set(gcf, 'UserData', data)
saveas(fH, fullfile(save_pth,'Mesh_difGratNoise_gamma_SNR_across_subjects_BinarizedDenoised_4'), 'fig');

%% SNR Mesh per stimulus type Broadband
threshold = 0;%3;
% gaussian weight for each stimuli
fH = figure(1005); clf, set(fH, 'name', 'Broadband weight' )

for c = 1:9
    subplot(3,3,c)
    data_to_plot = squeeze(nanmean(snr_w_pwr_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]* 2)
end

data = {};
data.snr_w_gauss   = snr_w_pwr;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};


hgexport(fH, fullfile(save_pth,'Mesh_stimulus_type_broadband_SNR_across_subjects_BinarizedDenoised_4'));
saveas(fH, fullfile(save_pth,'Mesh_stimulus_type_broadband_SNR_across_subjects_BinarizedDenoised_4'), 'fig');

%% SNR Mesh
threshold = 0;%3;
% gaussian weight for each stimuli
fH = figure(455); clf, set(fH, 'name', 'Broadband weight')
%for c = 1:12
for c = 10:13
    subplot(2,2,c-9)
    data_to_plot = squeeze(nanmean(snr_w_pwr_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]* 2)
end

data = {};
data.snr_w_gauss   = snr_w_pwr;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};

hgexport(fH, fullfile(save_pth,'Mesh_difGratNoise_broadband_SNR_across_subjects_BinarizedDenoised_4'));
saveas(fH, fullfile(save_pth,'Mesh_difGratNoise_broadband_SNR_across_subjects_BinarizedDenoised_4'), 'fig');

