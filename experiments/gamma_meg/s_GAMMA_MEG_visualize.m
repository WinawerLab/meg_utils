% s_GAMMA_MEG_visualize

% Visualizes summary statistics computed in s_GAMMA_MEG_analysis.m

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';


fs                            = 1000;
intertrial_trigger_num        = 10;
<<<<<<< HEAD
which_data_to_visualize       = 6;
=======
which_data_to_visualize       = 9;%4:6;
>>>>>>> 00065d3efb95d2d6fe31a5ee9eaf801c5980be8f
save_images                   = true;
using_denoised_data           = true;

% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);



%% Loop over data sets
for subject_num = which_data_to_visualize
    
    condition_names               = {   ...
        'White Noise' ...
        'Binary White Noise' ...
        'Pink Noise' ...
        'Brown Noise' ...
        'Gratings(0.36 cpd)' ...
        'Gratings(0.73 cpd)' ...
        'Gratings(1.45 cpd)' ...
        'Gratings(2.90 cpd)' ...
        'Plaid'...
        'Blank'};
    
    if subject_num >= 9
        condition_names{3} = 'Binary Pink Noise';
        condition_names{4} = 'Binary Brown Noise';
    end
    
    %% Define paths and load data
    load_pth    = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    
    if using_denoised_data
        save_pth = fullfile(project_pth,'/Images',subj_pths{subject_num}, 'denoised');
        d        =  dir(fullfile(load_pth, '*_denoisedData_*boot*'));
        badChannels = [];
    else
        save_pth = fullfile(project_pth,'/Images',subj_pths{subject_num});
        d        =  dir(fullfile(load_pth, '*boot*'));
        badChannels = zeros(1,157);
    end
    
    results         = load(fullfile(load_pth, d(1).name));
    
    num_conditions = size(results.out_exp,2);
    num_channels   = size(results.out_exp,1);
    
    if isempty(badChannels)
        denoisedData   = load(fullfile(project_pth, subj_pths{subject_num}, 'processed',sprintf('s0%d_denoisedData.mat',subject_num+1)));
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
    
    snr_out_exp = zeros(num_channels, num_contrasts);
    snr_w_gauss = zeros(num_channels, num_contrasts);
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
    threshold = 2;%3;
    % gaussian weight for each stimuli
    fH = figure(1005); clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);  set(fH, 'name', 'Gaussian SNR' )
    %for c = 1:12
    for c = 1:9
        subplot(3,3,c)
        data_to_plot = snr_w_gauss(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot,~badChannels,0), contrastnames{c});
        set(gca, 'CLim', [-1 1]* 3)
    end

    
    if save_images
        if ~exist(save_pth, 'dir'), mkdir(save_pth); end
        hgexport(fH, fullfile(save_pth,'Per_Condition_Gamma_SNR_thresh1.eps')); 
    end
    
    %% SNR Mesh for Broadband
    scrsz = get(0,'ScreenSize');
    threshold = 2;%3;
    % gaussian weight for each stimuli
    fH = figure(455); clf, set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); set(fH, 'name', 'Broadband SNR')
    %for c = 1:12
    for c = 1:9
        subplot(3,3,c)
        data_to_plot = snr_w_pwr(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(to157chan(data_to_plot,~badChannels,0), contrastnames{c});
        set(gca, 'CLim', [-1 1]* 3)
    end
    if save_images;  hgexport(fH, fullfile(save_pth,'Per_Condition_BB_SNR_thresh2')); end
   
    

    %% Noise (Gamma SNR - baseline)
    
    threshold = 2;
    fH = figure(2); clf, set(fH, 'name', 'Gamma weight')
    % noise gamma
    data_to_plot = snr_w_gauss(:,12);
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    
    [fH,ch] = megPlotMap(to157chan(data_to_plot',~badChannels,0),[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    title(contrastnames{12})
    
    if save_images; hgexport(fH, fullfile(save_pth,'figure_gammapower_noise_m_gratings_thresh2.eps')); end
    
    
    %% Gratings (gamma SNR - Baseline)
    fH = figure(2); clf, set(fH, 'name', 'Gamma weight')
    % Gratings gamma
    data_to_plot = snr_w_gauss(:,2)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    
    [fH,ch] = megPlotMap(to157chan(data_to_plot,~badChannels,0),[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    title(contrastnames{2})
    
    
    if save_images;  hgexport(fH, fullfile(project_pth,'figure_gammapower_gratings.eps')); end
    %% Gratings (Broadband SNR - Baseline)
    
    fH = figure(3); clf, set(fH, 'name', 'Broadband weight')
    %     for c = [2 1]
    %         subplot(1,2,c)
    threshold = 2;
    data_to_plot = to157chan(snr_w_pwr(:,12)',~badChannels,0);
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    %         ft_plotOnMesh(data_to_plot, contrastnames{c});
    %         set(gca, 'CLim', [-1 1]* 10)
    %     end
    
    [fH,ch] = megPlotMap(data_to_plot,[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{12})
    
    
    
    if save_images; hgexport(fH, fullfile(save_pth,'figure_bbpower_noise_m_gratings_thresh2.eps')); end
    
    
    
    
    
    
    
    %% Noise (Broadband SNR - Baseline)
    
    %     fH = figure(4); clf, set(fH, 'name', 'Broadband weight')
    %     for c = [2 1]
    %         subplot(1,2,c)
    data_to_plot = snr_w_pwr(:,1)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    %         ft_plotOnMesh(data_to_plot, contrastnames{c});
    %         set(gca, 'CLim', [-1 1]* 10)
    %     end
    
    [fH,ch] = megPlotMap(to157chan(data_to_plot,~badChannels,0),[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{1})
    
    hgexport(fH, fullfile(project_pth,'figure_bbpower_noise_s1.eps'));
    
    
    %% Meshes of Gaussian/Broadband Weight
    
    % TODO: threshold maps by significance: w_gauss_mn./w_gauss_sd>2
    
    
    fH = figure(998); clf, set(fH, 'name', 'Gaussian weight (mediann-baseline) before denoising')
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan((w_gauss_mn(:,cond)' - w_gauss_mn(:,num_conditions)'),~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .2])
    end
    
    if save_images
        figurewrite([],fH,[], fullfile(save_pth, 'Mesh_Gaussian_weight'));
    end
    
    fH = figure(999); clf, set(fH, 'name', 'Broadband weight')
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .1])
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Broadband.eps'));
        
    end
    
    
    
    %axis tight
    
    %% Mesh visualization of model fits
    
    % TODO: threshold maps by significance: w_gauss_mn./w_gauss_sd>2
    
    
    fH = figure(998); clf, set(fH, 'name', 'Gaussian weight'); set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan(w_gauss_mn(:,cond)' - w_gauss_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .1])
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Per_Condition_Gamma_weight(median-baseline).eps'));
    end
    
    fH = figure(999); clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(to157chan(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)',~badChannels,0), condition_names{cond});
        set(gca, 'CLim', [0 .1])
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Per_Condition_BB_weight(median-baseline).eps'));
    end
    
    %     fH = figure(1000); clf
    %     for cond = 1:9
    %         subplot(3,3,cond)
    %         ft_plotOnMesh(w_pwr_mn(:,cond)', condition_names{cond});
    %         set(gca, 'CLim', [0 2])
    %     end
    
    
    scrsz = get(0,'ScreenSize');
    fH = figure(1000); clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    subplot(2,2,1)
    ft_plotOnMesh(to157chan((w_gauss_mn * [0 0 0 0 1 1 1 1 0 -4]')',~badChannels,0), 'Gamma power, All gratings minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,2)
    ft_plotOnMesh(to157chan((w_gauss_mn * [1 1 1 1 0 0 0 0 0 -4]')',~badChannels,0), 'Gamma power, All noise minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,3)
    ft_plotOnMesh(to157chan((w_pwr_mn * [0 0 0 0 1 1 1 1 0 -4]')',~badChannels,0), 'Broadband power, All gratings minus baseline');
    set(gca, 'CLim', .25*[-1 1])
    
    subplot(2,2,4)
    ft_plotOnMesh(to157chan((w_pwr_mn * [1 1 1 1 0 0 0 0 0 -4]')',~badChannels,0), 'Broadband power, All noise minus baseline');
    set(gca, 'CLim', .25*[-1 1])
    
    if save_images
        hgexport(1000, fullfile(save_pth, 'Mesh_G_BB_Noise_or_Gratings_M_Baseline.eps'));
    end
    
    fH = figure(1001); clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    subplot(2,2,1)
    ft_plotOnMesh(to157chan((w_gauss_mn * [-1 -1 -1 -1 1 1 1 1 0 0]')',~badChannels,0), 'Gamma power, All gratings minus all noise');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,2)
    ft_plotOnMesh(to157chan((w_pwr_mn * [-1 -1 -1 -1 1 1 1 1 0 0]')',~badChannels,0), 'Broadband power, All gratings minus all noise');
    set(gca, 'CLim', .5*[-1 1])
    
    
    subplot(2,2,3)
    ft_plotOnMesh(to157chan((w_gauss_mn * [1 1 1 1 1 1 1 1 1 -9]')',~badChannels,0), 'Gamma power, All stimuli minus baseline');
    set(gca, 'CLim', .5 *[-1 1])
    
    subplot(2,2,4)
    ft_plotOnMesh(to157chan((w_pwr_mn * [1 1 1 1 1 1 1 1 1 -9]')',~badChannels,0), 'Broadband, All stimuli minus baseline');
    set(gca, 'CLim', .5 * [-1 1])
    
    if save_images
        hgexport(1001, fullfile(save_pth, 'Mesh_G_BB_Gratings_M_noise_All_M_Baseline.eps'));
    end
    
 end