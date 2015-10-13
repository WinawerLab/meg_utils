% s_GAMMA_MEG_visualize

% Visualizes summary statistics computed in s_GAMMA_MEG_analysis.m


fs                            = 1000;
intertrial_trigger_num        = 10;
which_session_to_visualize    = 5:6; %[7:12,14:16];
save_images                   = true;
using_denoised_data           = false;
suffix                        = 'localregression_multi_100';

if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end




%% Loop over data sets
for session_num = which_session_to_visualize
    
   condition_names = gamma_get_condition_names(session_num);

    path_to_data = meg_gamma_get_path(session_num);
    %% Define paths and load data
    load_pth    = fullfile(path_to_data, 'processed');
    
    if using_denoised_data
        save_pth = fullfile(project_pth,'/Images',subj_pths{session_num-1}, 'denoised');
        d        =  dir(fullfile(load_pth, '*_denoisedData_*boot*'));tmp
        badChannels = [];
    else
        save_pth = fullfile(project_pth,'/Images',subj_pths{session_num-1});
        d        =  dir(fullfile(load_pth, sprintf('*%s*', suffix)));
        badChannels = zeros(1,157);
    end
    
    results         = load(fullfile(load_pth, d(1).name));
    
    num_conditions = size(results.w_pwr,2);
    num_channels   = size(results.w_pwr,1);
    
    if isempty(badChannels)
        denoisedData   = load(fullfile(project_pth, subj_pths{session_num-1}, 'processed',sprintf('s0%d_denoisedData.mat',session_num)));
        badChannels    = denoisedData.bad_channels;
    end
    
    
    %% Calculating SNR contrasts
    if results.nboot > 1,
        summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
    else
        summary_stat = @(x) nanmean(x,3);
    end
    
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
        sprintf('%s - %s', condition_names{1}, condition_names{10})...
        sprintf('%s - %s', condition_names{2}, condition_names{10})...
        sprintf('%s - %s', condition_names{3}, condition_names{10})...
        sprintf('%s - %s', condition_names{4}, condition_names{10})...
        sprintf('%s - %s', condition_names{5}, condition_names{10})...
        sprintf('%s - %s', condition_names{6}, condition_names{10})...
        sprintf('%s - %s', condition_names{7}, condition_names{10})...
        sprintf('%s - %s', condition_names{8}, condition_names{10})...
        sprintf('%s - %s', condition_names{9}, condition_names{10})...        
        'noise - baseline'...
        'gratings - baseline'...
        'noise - gratings'...
        'gratings - noise'...
        };
    
    
    % ensure each condition is weighted proportionatly in each contrast
    contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
    
    num_contrasts = size(contrasts,1);
    
    % compute means
    
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
    
    % threshold (replace SNR values < 2 or > 20 with 0)
    
    %% SNR Mesh for Gaussian bump
    scrsz = get(0,'ScreenSize');
    threshold = 0;%3;
    
    % gaussian weight for each stimuli
    fH = figure; clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);  set(fH, 'name', 'Gaussian SNR' )
    plot_range = [-1 1] * ceil(max(max(abs(snr_w_gauss(:,1:9)))));
    for c = 1:9
        subplot(3,3,c)
        data_to_plot = snr_w_gauss(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot,~badChannels,0), contrastnames{c});
        set(gca, 'CLim', plot_range)
        colormap parula
    end
    
    if save_images
        if ~exist(save_pth, 'dir'), mkdir(save_pth); end
        hgexport(fH, fullfile(save_pth,sprintf('Per_Condition_Gamma_SNR_local_%s.eps',suffix))); 
        close(fH)
    end
    
    %% SNR Mesh for Broadband
    scrsz = get(0,'ScreenSize');
    threshold = 0;%3;
    fH = figure; clf, set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); set(fH, 'name', 'Broadband SNR')
    plot_range = [-1 1] * ceil(max(max(abs(snr_w_pwr(:,1:9)))));

    for c = 1:9
        subplot(3,3,c)
        data_to_plot = snr_w_pwr(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot,~badChannels,0), contrastnames{c});
        set(gca, 'CLim', plot_range)
        colormap parula
    end
    if save_images;  hgexport(fH, fullfile(save_pth,sprintf('Per_Condition_BB_SNR_%s',suffix))); close(fH);end
   
    

    %     %% Noise (Gamma SNR - baseline)
    %     this_contrast = 12;
    %
    %     threshold = 0;
    %     fH = figure;  set(fH, 'name', 'Gamma weight')
    %     % noise gamma
    %     data_to_plot = snr_w_gauss(:,this_contrast);
    %     data_to_plot(abs(data_to_plot) < threshold) = 0;
    %     plot_range = [-1 1]*max(abs(data_to_plot));
    %     [~,ch] = megPlotMap(to157chan(data_to_plot',~badChannels,0),plot_range,gcf,jmaColors('coolhotcortex'));
    %
    %     makeprettyaxes(gca,9,9);
    %     set(ch,'ytick',-10:10:10);
    %     title(contrastnames{this_contrast})
    %
    %     if save_images; hgexport(fH, fullfile(save_pth,sprintf('figure_gammapower_noise_m_gratings_%s.eps',suffix))); end
    %
    %
    %     %% Gratings (gamma SNR - Baseline)
    %     this_contrast = 13;
    %     fH = figure; clf, set(fH, 'name', 'Gamma weight')
    %     % Gratings gamma
    %     data_to_plot = snr_w_gauss(:,this_contrast)';
    %     data_to_plot(abs(data_to_plot) < threshold) = 0;
    %     plot_range = [-1 1]*max(abs(data_to_plot));
    %
    %     [fH,ch] = megPlotMap(to157chan(data_to_plot,~badChannels,0),plot_range,gcf,jmaColors('coolhotcortex'));
    %
    %     makeprettyaxes(gca,9,9);
    %     set(ch,'ytick',-10:10:10);
    %     title(contrastnames{this_contrast})
    %
    %
    %     if save_images;  hgexport(fH, fullfile(project_pth,sprintf('figure_gammapower_gratings_%s.eps',suffix))); end
    %     %% Gratings - Baseline (Broadband SNR)
    %     this_contrast = 11;
    %
    %     fH = figure; clf, set(fH, 'name', 'Broadband weight')
    %     threshold = 0;
    %     data_to_plot = to157chan(snr_w_pwr(:,this_contrast)',~badChannels,0);
    %     data_to_plot(abs(data_to_plot) < threshold) = 0;
    %     plot_range = [-1 1]*max(abs(data_to_plot));
    %
    %     [fH,ch] = megPlotMap(data_to_plot,plot_range,gcf,jmaColors('coolhotcortex'));
    %
    %     makeprettyaxes(gca,9,9);
    %     title(contrastnames{this_contrast})
    %
    %     if save_images; hgexport(fH, fullfile(save_pth,sprintf('figure_bbpower_noise_m_gratings_thresh2_%s.eps',suffix))); end
    %
    %
    
    %% Meshes of Gaussian/Broadband Weight
    
    % TODO: threshold maps by significance: w_gauss_mn./w_gauss_sd>2
    
    
%     fH = figure(998); clf, set(fH, 'name', 'Gaussian weight (median-baseline) before denoising')
%     for cond = 1:9
%         subplot(3,3,cond)
%         ft_plotOnMesh(to157chan((w_gauss_mn(:,cond)' - w_gauss_mn(:,num_conditions)'),~badChannels,0), condition_names{cond});
%         set(gca, 'CLim', [0 .2])
%         colormap parula
% 
%     end
%     
%     if save_images
%         hgexport(fH,fullfile(save_pth, 'Mesh_Gaussian_weight_per_condition_local_regression1.eps'));
%     end
%     
%     fH = figure(999); clf, set(fH, 'name', 'Broadband weight')
%     for cond = 1:9
%         subplot(3,3,cond)
%         ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
%         set(gca, 'CLim', [0 .2])
%         colormap parula
% 
%     end
%     
%     if save_images
%         hgexport(fH, fullfile(save_pth, sprintf('Mesh_Broadband_per_condition_%s.eps',suffix)));
%         
%     end
    
    
    
    %axis tight
    
    %% Mesh visualization of model fits
    
    % TODO: threshold maps by significance: w_gauss_mn./w_gauss_sd>2
    
    
    fH = figure; clf, set(fH, 'name', 'Gaussian weight'); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan(w_gauss_mn(:,cond)' - w_gauss_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .2])
        colormap parula

    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, sprintf('Per_Condition_Gamma_weight(median-baseline)_%s.eps',suffix)));
        close(fH)
    end
    
    fH = figure; clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .2])
        colormap parula

    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, sprintf('Per_Condition_BB_weight(median-baseline)_%s.eps',suffix)));
        close(fH)
    end
    
   
    
 end