%% s_GAMMA_MEG_visualize_natural
% visualize processed data from s_GAMMA_MEG_analysis.m
% TODO: loop across data sets

%% Parameters

% paths
if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end
project_pth    = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth       = '*_Gamma_*subj*';
% Find subject path
d = dir(fullfile(project_pth, data_pth));
%   restrict to directories
subj_pths = struct2cell(d);

% data parameters
fs                            = 1000;
intertrial_trigger_num        = 14;
which_session_to_visualize    = 18; %[7:12,14:16];
save_images                   = false;
using_denoised_data           = true;
suffix                        = 'localregression_multi_100_boots';
conditions                    = gamma_get_condition_names(which_session_to_visualize);
 %% loop over sessions
for session_num = which_session_to_visualize
    %% load data
    condition_names = gamma_get_condition_names(session_num);
    path_to_data = meg_gamma_get_path(session_num);
    [~, subject_folder] = fileparts(path_to_data);

    load_pth    = fullfile(path_to_data, 'processed');
    
    if using_denoised_data
        save_pth = fullfile(project_pth,subject_folder,'figs','denoised');
        d        =  dir(fullfile(load_pth, '*_denoisedData_*boot*'));
        badChannels = zeros(1,157);
    else
        save_pth = fullfile(project_pth,subject_folder,'figs');
        d        =  dir(fullfile(load_pth, sprintf('*%s*', suffix)));
        badChannels = zeros(1,157);
    end
    
    % check this because in the previous script d(2) pathed to the summary
    % stats
    results         = load(fullfile(load_pth, d(1).name)); % f x cond x chan
    
    %% calculate contrasts
    
    num_conditions = size(results.w_pwr,2);
    num_channels   = size(results.w_pwr,1);
    
    
    % take the mean across channels
    if results.nboot > 1,
        summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
    else
        summary_stat = @(x) nanmean(x,3);
    end
    
    % difference identity vectors
    % each condition - baseline
    contrasts = eye(num_conditions,num_conditions); %eye cond x cond
    contrasts(:,13) = -1; % make baseline column negative to subtract from each cond
    
    contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
    
    num_contrasts = size(contrasts,1);
    
    w_gauss_mn = nanmean(results.w_gauss,3);
    w_pwr_mn   = nanmean(results.w_pwr,3);
    
    % compute SNR
    snr_fit_bl  = zeros(num_channels, num_contrasts);
    snr_gauss_f = zeros(num_channels, num_contrasts);
    
    tmp_data = permute(results.w_pwr, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, results.nboot);
    snr_w_pwr = summary_stat(tmp)';
    
    tmp_data = permute(results.w_gauss, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, results.nboot);
    snr_w_gauss  = summary_stat(tmp)';
    
    %% mesh plots
    % SNR Gauss
    
    scrsz = get(0,'ScreenSize');
    threshold = 0;%3;
    
    % gaussian weight for each stimuli
    fH = figure; clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)]);  set(fH, 'name', 'Gaussian SNR' )
    plot_range = [-1 1] * (max(max(abs(snr_w_gauss(:,1:length(conditions))))));
    for c = 1:12
        subplot(4,3,c)
        data_to_plot = snr_w_gauss(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot,badChannels,0), conditions{c});
        set(gca, 'CLim', plot_range)
        colormap parula
    end
    
    if save_images
        if ~exist(save_pth, 'dir'), mkdir(save_pth); end
        hgexport(fH, fullfile(save_pth,sprintf('%s_Per_Condition_Gamma_SNR_local_%s.eps',date,suffix)));
        close(fH)
    end
    
    %%
    %SNR Power
    
    scrsz = get(0,'ScreenSize');
    threshold = 0;%3;
    fH = figure; clf, set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)]); set(fH, 'name', 'Broadband SNR')
    plot_range = [-1 1] * (max(max(abs(snr_w_pwr(:,1:length(conditions))))));
    
    for c = 1:12
        subplot(4,3,c)
        data_to_plot = snr_w_pwr(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot',~badChannels,0), conditions{c});
        set(gca, 'CLim', plot_range)
        colormap parula
    end
    if save_images;  
        hgexport(fH, fullfile(save_pth,sprintf('%s_Per_Condition_BB_SNR_%s.eps',date, suffix))); 
        close(fH);
    end
    
end
