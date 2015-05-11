% s_GAMMA_MEG_visualize

% Visualizes summary statistics computed in s_GAMMA_MEG_analysis.m

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';

data_channels                 = 1:157;
num_channels                  = length(data_channels);
fs                            = 1000;
intertrial_trigger_num        = 10;
which_data_to_visualize       = 5;
save_images                   = false;

meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

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

%% Loop over data sets
for subject_num = which_data_to_visualize
    %% Load Data
    load_pth    = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    d           =  dir(fullfile(load_pth, '*boot*'));
    res         = load(fullfile(load_pth, d(1).name));
    
    num_conditions = size(res.out_exp,2);
    
    %% Calculating SNR contrasts
    if res.nboot > 1,
        summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
    else
        summary_stat = @(x) nanmean(x,3);
    end
    
    contrasts = [...
        1 1 1 1 0 0 0 0 0 -4; ...    % noise - baseline
        0 0 0 0 1 1 1 1 0 -4; ...    % gratings - baseline
        1 1 1 1 -1 -1 -1 -1 0 0; ... % noise - gratings
        -1 -1 -1 -1 1 1 1 1 0 0; ... % gratings - noise
        1 0 0 0 0 0 0 0 0 -1; ...    % white noise - baseline
        0 1 0 0 0 0 0 0 0 -1; ...    % binarized white noise - baseline
        0 0 1 0 0 0 0 0 0 -1; ...    % pink noise - baseline
        0 0 0 1 0 0 0 0 0 -1; ...    % brown noise - baseline
        0 0 0 0 1 0 0 0 0 -1; ...    % 0.36cpd gratings - baseline
        0 0 0 0 0 1 0 0 0 -1; ...    % 0.73cpd gratings - baseline
        0 0 0 0 0 0 1 0 0 -1; ...    % 1.46cpd gratings - baseline
        0 0 0 0 0 0 0 1 0 -1; ...    % 2.90cpd gratings - baseline
        0 0 0 0 0 0 0 0 1 -1];       % plaid - baseline
    
    contrastnames = {
        'noise - baseline'...
        'gratings - baseline'...
        'noise - gratings'...
        'gratings - noise'...
        'white noise - baseline'...
        'binwn - baseline'...
        'pink noise - baseline'...
        'brown noise - baseline'...
        '0.36cpd gratings - baseline'...
        '0.73cpd gratings - baseline'...
        '1.46cpd gratings - baseline'...
        '2.90cpd gratings - baseline'...
        'plaid - baseline'};
    
    
    % ensure each condition is weighted proportionatly in each contrast
    contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
    
    num_contrasts = size(contrasts,1);
    
    % compute means
    
    w_gauss_mn = nanmean(res.w_gauss,3);
    w_pwr_mn   = nanmean(res.w_pwr,3);
    
    % compute SNR
    
    snr_out_exp = zeros(num_channels, num_contrasts);
    snr_w_gauss = zeros(num_channels, num_contrasts);
    snr_gauss_f = zeros(num_channels, num_contrasts);
    
    tmp_data = permute(res.w_pwr, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, res.nboot);
    snr_w_pwr = summary_stat(tmp)';
    
    tmp_data = permute(res.w_gauss, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, num_channels, res.nboot);
    snr_w_gauss  = summary_stat(tmp)';
    
    % threshold (replace SNR values < 2 or > 20 with 0)
    
    %% SNR Mesh
    threshold = 0;%3;
    % gaussian weight for each stimuli
    fH = figure(998); clf, set(fH, 'name', 'Gaussian weight')
    for c = 1:12
        subplot(3,4,c)
        data_to_plot = snr_w_gauss(:,c)';
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(data_to_plot, contrastnames{c});
        set(gca, 'CLim', [-1 1]* 10)
    end
    
    %% Noise (Gamma SNR - baseline)
    
    threshold = 1;
    %     fH = figure(1); clf, set(fH, 'name', 'Gamma weight')
    % noise gamma
    data_to_plot = snr_w_gauss(:,1)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    
    [fH,ch] = megPlotMap(data_to_plot,[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{1})
    
    hgexport(fH, fullfile(project_pth,'figure_gammapower_noise_s1.eps'));
    
    
    %% Gratings (gamma SNR - Baseline)
    %         fH = figure(2); clf, set(fH, 'name', 'Gamma weight')
    % Gratings gamma
    data_to_plot = snr_w_gauss(:,2)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    
    [fH,ch] = megPlotMap(data_to_plot,[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{2})
    
    
    hgexport(fH, fullfile(project_pth,'figure_gammapower_gratings_s1.eps'));
    %% Gratings (Broadband SNR - Baseline)
    
    %     fH = figure(3); clf, set(fH, 'name', 'Broadband weight')
    %     for c = [2 1]
    %         subplot(1,2,c)
    data_to_plot = snr_w_pwr(:,2)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    %         ft_plotOnMesh(data_to_plot, contrastnames{c});
    %         set(gca, 'CLim', [-1 1]* 10)
    %     end
    
    [fH,ch] = megPlotMap(data_to_plot,[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{2})
    
    
    
    hgexport(fH, fullfile(project_pth,'figure_bbpower_gratings_s1.eps'));
    
    
    
    
    
    
    
    %% Noise (Broadband SNR - Baseline)
    
    %     fH = figure(4); clf, set(fH, 'name', 'Broadband weight')
    %     for c = [2 1]
    %         subplot(1,2,c)
    data_to_plot = snr_w_pwr(:,1)';
    data_to_plot(abs(data_to_plot) < threshold) = 0;
    %         ft_plotOnMesh(data_to_plot, contrastnames{c});
    %         set(gca, 'CLim', [-1 1]* 10)
    %     end
    
    [fH,ch] = megPlotMap(data_to_plot,[-10,10],gcf,jmaColors('coolhotcortex'));
    
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:10:10);
    %     makeprettyaxes(ch,9,9);
    title(contrastnames{1})
    
    hgexport(fH, fullfile(project_pth,'figure_bbpower_noise_s1.eps'));
    
    %% Plot Gaussian fits
    
    % include spectral_data_boots to the HPC save script to produce
    % spectral_data_mean
    
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
            plot(f(f_sel),data_fit(f_sel),'b-','LineWidth',line_width)
            
            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                'XLim', [min(f_use4fit) max(f_use4fit)], ...
                'YLim', 10.^[0.3 1.5] )
            title(sprintf('Subject %d, Channel %d, %s', subject_num, chan, condition_names{cond}))
        end
        if save_images,
            
            hgexport(fH, fullfile(save_pth, sprintf('Spectra_Chan%03d.eps', chan)));
        else
            waitforbuttonpress
        end
    end
    
    %% Plot spectra
    fH = figure(2); clf,  set(gcf, 'Position', [100 100 1000 400], 'Color', 'w')
    
    colors = zeros(num_conditions,3);
    colors(1:4,:) = ([0 1 0]'*[.4 .6 .8 1])';
    colors(5:9,:) = ([1 0 0]'*[.4 .55 .7 .85 1])'; %hsv(5);
    colors(10,:)  = [0 0 0];
    yl            =  10.^([0.5 1.6]);
    xl            = [30 200];
    lw = 4;
    for chan = data_channels
        
        % Plot Data -----------------------------------------------------------
        subplot(1,3,1); cla; hold on
        for cond = 1:num_conditions
            plot(f(f_sel), squeeze(spectral_data_mean(f_sel,cond,chan)), 'Color', colors(cond,:), 'LineWidth', lw);
        end
        set(gca, 'YScale', 'log','XScale', 'log', 'YLim', yl, 'XLim', xl,  ...
            'Color', [1 1 1], 'XGrid', 'on', ...
            'XTick', [10 60 100 200]);
        title(sprintf('Data from channel %d, subject %d', chan, subject_num))
        
        % Plot Fits -----------------------------------------------------------
        subplot(1,3,2); cla; hold on
        for cond = 1:num_conditions
            plot(f, squeeze(10.^(fit_f2_mn(cond,:, chan))), 'Color', colors(cond,:), 'LineWidth', lw);
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
        else
            waitforbuttonpress
        end
        
    end
    %axis tight
    
    %% Meshes of Gaussian/Broadband Weight
    
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
    
    fH = figure(999); clf, set(fH, 'name', 'Broadband weight')
    for cond = 1:9
        subplot(3,3,cond)
        ft_plotOnMesh(w_pwr_mn(:,cond)' - w_pwr_mn(:,num_conditions)', condition_names{cond});
        set(gca, 'CLim', [-1 1] *.03)
    end
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Broadband.eps'));
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
    ft_plotOnMesh((w_gauss_mn * [0 0 0 0 1 1 1 1 0 -4]')', 'Gamma power, All gratings minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,2)
    ft_plotOnMesh((w_gauss_mn * [1 1 1 1 0 0 0 0 0 -4]')', 'Gamma power, All noise minus baseline');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,3)
    ft_plotOnMesh((w_gauss_mn * [-1 -1 -1 -1 1 1 1 1 0 0]')', 'Gamma power, All gratings minus all noise');
    set(gca, 'CLim', .5*[-1 1])
    
    subplot(2,2,4)
    ft_plotOnMesh((w_pwr_mn * [1 1 1 1 1 1 1 1 1 -9]')', 'Broadband, All stimuli minus baseline');
    set(gca, 'CLim', .2 * [-1 1])
    
    if save_images
        hgexport(fH, fullfile(save_pth, 'Mesh_Gamma_Gratings_M_Baseline.eps'));
        
    end
    
end