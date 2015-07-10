fH = figure; set(fH, 'Color', 'w');


% Visualizes summary statistics computed in s_GAMMA_MEG_analysis.m

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';


fs                            = 1000;
intertrial_trigger_num        = 10;
which_data_to_visualize       = 4:6;
save_images                   = true;

meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

get_snr  = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
get_mean = @(x) nanmean(x,3);
get_std  = @(x) nanstd(x,[],3);

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
num_conditions = 10;
nboot = 100;


%% Loop over data sets
for subject_num = which_data_to_visualize
    
    load_pth        = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    datasets        = dir(fullfile(load_pth, '*boot*'));
    before          = load(fullfile(load_pth, datasets(1).name));
    save_pth        = fullfile(project_pth,'/Images',subj_pths{subject_num}, 'denoised');
    if subject_num < 6
        after       = load(fullfile(load_pth, datasets(3).name));
    else
        after       = load(fullfile(load_pth, datasets(2).name));
    end
    
    denoisedData = dir(fullfile(load_pth, sprintf('s0%d_denoisedData*',subject_num+1)));
    denoisedData = load(fullfile(load_pth,denoisedData(1).name));
    
    % compute S, N, SNR Gamma
    tmp_data = permute(before.w_gauss(~denoisedData.bad_channels,:,:), [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, size(before.w_gauss(~denoisedData.bad_channels,:,:),1), nboot);
    signal_before_gamma = get_mean(tmp)';
    noise_before_gamma = get_std(tmp)';
    snr_before_gamma = get_snr(tmp)';
    
    tmp_data = permute(after.w_gauss, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, size(after.w_gauss,1), nboot);
    signal_after_gamma = get_mean(tmp)';
    noise_after_gamma = get_std(tmp)';
    snr_after_gamma = get_snr(tmp)';
    
    % compute S, N, SNR Broadband
    tmp_data = permute(before.w_pwr(~denoisedData.bad_channels,:,:), [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, size(before.w_pwr(~denoisedData.bad_channels,:,:),1), nboot);
    signal_before_bb = get_mean(tmp)';
    noise_before_bb = get_std(tmp)';
    snr_before_bb  = get_snr(tmp)';
    
    tmp_data = permute(after.w_pwr, [2 1 3]);
    tmp_data = reshape(tmp_data, num_conditions, []);
    tmp = contrasts*tmp_data;
    tmp = reshape(tmp, num_contrasts, size(after.w_pwr,1), nboot);
    signal_after_bb = get_mean(tmp)';
    noise_after_bb = get_std(tmp)';
    snr_after_bb  = get_snr(tmp)';

    
    %     % max across 3 conditions
    %     finalsnr = cat(2, snr_before_gamma, snr_after_gamma);
    %     % max across before and after
    %     finalsnr = max(finalsnr');
    %     % sort
    %     [~,idx] = sort(finalsnr,'descend');
    %     % find the top 10
    %     pcchan_gamma = false(size(denoisedData.results.noisepool));
    %     pcchan_gamma(idx(1:10))= 1;
    
    
    % signal and noise before denoising
    g_signal1 = signal_before_gamma;
    g_noise1  = noise_before_gamma;
    % signal and noise after denoising
    g_signal2 = signal_after_gamma;
    g_noise2  =  noise_after_gamma;
    % snr = signal/denoise
    g_snr1    = snr_before_gamma;
    g_snr2    = snr_after_gamma;
    
    
    % signal and noise before denoising
    bb_signal1 = signal_before_bb;
    bb_noise1  = noise_before_bb;
    % signal and noise after denoising
    bb_signal2 = signal_after_bb;
    bb_noise2  =  noise_after_bb;
    % snr = signal/denoise
    bb_snr1    = snr_before_bb;
    bb_snr2    = snr_after_bb;
    
    condColors = jet(9);
    
    % GAMMA plot each condition as a different color
    % signal
    figure(1); set(gcf,'color','w');
    subject_num2 = subject_num-3;
    subplot(3,3,subject_num2); cla; hold on;
    for nn = 1:9
        plot(g_signal1(:,nn),g_signal2(:,nn),'o','color',condColors(nn,:), 'MarkerSize',8, 'LineWidth',2);
    end
    axis square;
    axismax = max([g_signal1(:);g_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    title(sprintf('S%d : signal', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    
    % noise
    subplot(3,3,subject_num2+length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(g_noise1(:,nn),g_noise2(:,nn),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([g_noise1(:); g_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    axis square;
    title(sprintf('S%d : noise', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    % snr
    subplot(3,3,subject_num2+2*length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(g_snr1(:,nn),g_snr2(:,nn),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([g_snr1(:); g_snr2(:)])*1.2;
    xlim([0,axismax]);  set(gca,'FontSize',18); ylim([0,axismax]);  set(gca,'FontSize',18); line([0,axismax],[0,axismax],'color','k', 'LineWidth',2);
    set(gca,'FontSize',18)
    axis square;
    title(sprintf('S%d : SNR', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    drawnow;
    
    hgexport(gcf, fullfile(save_pth, 'S_N_SNR_scatterplot_Gamma.eps'));
    %% BROADBAND
    % GAMMA plot each condition as a different color
    % signal
    figure(2); set(gcf,'color','w');
    subplot(3,3,subject_num2); cla; hold on;
    for nn = 1:9
        plot(bb_signal1(:,nn),bb_signal2(:,nn),'o','color',condColors(nn,:), 'MarkerSize',8, 'LineWidth',2);
    end
    axis square;
    axismax = max([bb_signal1(:);bb_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    title(sprintf('S%d : signal', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    
    % noise
    subplot(3,3,subject_num2+length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(bb_noise1(:,nn),bb_noise2(:,nn),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([bb_noise1(:); bb_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    axis square;
    title(sprintf('S%d : noise', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    % snr
    subplot(3,3,subject_num2+2*length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(bb_snr1(:,nn),bb_snr2(:,nn),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([bb_snr1(:); bb_snr2(:)])*1.2;
    xlim([0,axismax]);  set(gca,'FontSize',18); ylim([0,axismax]);  set(gca,'FontSize',18); line([0,axismax],[0,axismax],'color','k', 'LineWidth',2);
    set(gca,'FontSize',18)
    axis square;
    title(sprintf('S%d : SNR', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    drawnow;
    
    hgexport(gcf, fullfile(save_pth, 'S_N_SNR_scatterplot_Broadband.eps'));
end