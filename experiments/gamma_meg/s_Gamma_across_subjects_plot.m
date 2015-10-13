%% GAMMA MEG make across subjects headplot

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';

% Subjects in the binarized noise experiment are: [10:12,14:16]
% Subjects in the non-binarized noise experiment are: [5:9]


fs                            = 1000;
intertrial_trigger_num        = 10;
session_nums                  = [5:9]; %[10:12,14:16];
save_images                   = false;
save_pth                      = fullfile(project_pth, 'Images');
use_denoised_data             = false;
if use_denoised_data;
     delete_me = 4; 
else delete_me = 1; 
end;


suffix                         = 'NotBinarizedDenoised_localregression';

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

scrsz = get(0,'ScreenSize');


for session = session_nums
    

    
    %% Load Data
    load_pth                = fullfile(project_pth, subj_pths{session -1}, 'processed');
    if use_denoised_data; d =  dir(fullfile(load_pth, '*denoisedData_bootstrapped100_2*'));  % Check whether d(1) or d(2) is used, for positive or both sides modelfit
    else                  d =  dir(fullfile(load_pth, '*_localregression*100.mat')); end
        
    AllResults{session-delete_me} = load(fullfile(load_pth, d(1).name));
    w_gauss                 = AllResults{session -delete_me}.w_gauss;
    w_pwr                   = AllResults{session -delete_me}.w_pwr;
    
    % Load denoisedData to get bad channels.
    if use_denoised_data
        d                       =  dir(fullfile(load_pth, '*denoisedData.mat'));
        denoisedData            = load(fullfile(load_pth, d(1).name));
        bad_channels{session-delete_me} = denoisedData.bad_channels;
    else bad_channels{session-1} = zeros(1,157);
    end
    
    
    
    num_conditions          = size(AllResults{session-delete_me}.w_gauss,2);
    num_channels            = size(AllResults{session-delete_me}.w_gauss,1);
    
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
    tmp = reshape(tmp, num_contrasts, num_channels, AllResults{session-delete_me}.nboot);
    snr_w_pwr{session-delete_me} = summary_stat(tmp)';
    
    tmp_data = permute(w_gauss, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, AllResults{session-delete_me}.nboot);
    snr_w_gauss{session-delete_me}  = summary_stat(tmp)';
    
end

snr_w_gauss_157 = [];
snr_w_pwr_157   = [];
for ii = 1:num_contrasts
    for jj = session_nums-delete_me
        snr_w_gauss_157{ii,jj} = to157chan(snr_w_gauss{jj}(:,ii)',~bad_channels{jj}, 'nans');
        snr_w_pwr_157{ii,jj} = to157chan(snr_w_pwr{jj}(:,ii)',~bad_channels{jj}, 'nans');      
    end
end

snr_w_gauss_157 = reshape(catcell(1,snr_w_gauss_157),[num_contrasts,numel(session_nums),157]);
snr_w_pwr_157   = reshape(catcell(1,snr_w_pwr_157),[num_contrasts,numel(session_nums),157]);


%% SNR Mesh per stimulus type Gamma
threshold = 0;%3;
% gaussian weight for each stimuli
fH = figure(1); clf, set(fH, 'name', 'Gaussian weight' ); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

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

hgexport(fH, fullfile(save_pth,sprintf('Mesh_stimulus_type_gamma_SNR_across_subjects_%s',suffix)));
set(gcf, 'UserData', data);
saveas(fH, fullfile(save_pth,sprintf('Mesh_stimulus_type_gamma_SNR_across_subjects_%s',suffix)), 'fig');

%% SNR Mesh difference Noise and Gratings
threshold = 0;%3;
% gaussian weight for across stimuli
fH = figure(2); clf, set(fH, 'name', 'Gaussian weight'); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
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

hgexport(fH, fullfile(save_pth,sprintf('Mesh_difGratNoise_gamma_SNR_across_subjects_%s',suffix)));
set(gcf, 'UserData', data)
saveas(fH, fullfile(save_pth,sprintf('Mesh_difGratNoise_gamma_SNR_across_subjects_%s',suffix)), 'fig');

%% SNR Mesh per stimulus type Broadband
threshold = 0;%3;
% broadband weight for each stimuli
fH = figure(3); clf, set(fH, 'name', 'Broadband weight' ); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

for c = 1:9
    subplot(3,3,c)
    data_to_plot = squeeze(nanmean(snr_w_pwr_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]* 1.5)
end

data = {};
data.snr_w_gauss   = snr_w_pwr;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};


hgexport(fH, fullfile(save_pth,sprintf('Mesh_stimulus_type_broadband_SNR_across_subjects_%s',suffix)));
saveas(fH, fullfile(save_pth,sprintf('Mesh_stimulus_type_broadband_SNR_across_subjects_%s',suffix)), 'fig');

%% SNR Mesh
threshold = 0;%3;
% broadband weight for across stimuli
fH = figure(4); clf, set(fH, 'name', 'Broadband weight'); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
%for c = 1:12
for c = 10:13
    subplot(2,2,c-9)
    data_to_plot = squeeze(nanmean(snr_w_pwr_157(c,:,:),2))';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    ft_plotOnMesh(data_to_plot, contrastnames{c});
    set(gca, 'CLim', [-1 1]* 1.5)
end

data = {};
data.snr_w_gauss   = snr_w_pwr;
data.data_to_plot  = data_to_plot;
data.contrastnames = {contrastnames{c}};

hgexport(fH, fullfile(save_pth,sprintf('Mesh_difGratNoise_broadband_SNR_across_subjects_%s', suffix)));
saveas(fH, fullfile(save_pth,sprintf('Mesh_difGratNoise_broadband_SNR_across_subjects_%s',suffix)), 'fig');

