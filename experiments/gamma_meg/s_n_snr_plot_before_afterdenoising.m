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
    
    load_pth    = fullfile(project_pth, subj_pths{subject_num}, 'processed');
    datasets        =  dir(fullfile(load_pth, '*boot*'));
    before         = load(fullfile(load_pth, datasets(1).name));
    if subject_num < 6
        after         = load(fullfile(load_pth, datasets(3).name));
    else 
         after         = load(fullfile(load_pth, datasets(2).name));
    end
         
    denoisedData = dir(fullfile(load_pth, sprintf('s0%d_denoisedData*',subject_num+1)));
    denoisedData = load(fullfile(load_pth,denoisedData(1).name));
    
    signal_before_gamma = nanmedian(before.w_gauss(~denoisedData.bad_channels,:,:),3);
    noise_before_gamma = nanstd(before.w_gauss(~denoisedData.bad_channels,:,:),[],3);
    
    signal_after_gamma = nanmedian(after.w_gauss,3);
    noise_after_gamma = nanstd(after.w_gauss,[],3);
    
    snr_before_gamma = signal_before_gamma ./ noise_before_gamma;
    snr_after_gamma = signal_after_gamma ./ noise_after_gamma;
    
    % max across 3 conditions 
    finalsnr = cat(2, snr_before_gamma, snr_after_gamma);
    % max across before and after
    finalsnr = max(finalsnr');
    % sort
    [~,idx] = sort(finalsnr,'descend');
    % find the top 10
    pcchan_gamma = false(size(denoisedData.results.noisepool));
    pcchan_gamma(idx(1:10))= 1;
    
    
    % signal and noise before denoising
    g_signal1 = signal_before_gamma(pcchan_gamma,:);
    g_noise1  = noise_before_gamma(pcchan_gamma,:);
    % signal and noise after denoising 
    g_signal2 = signal_after_gamma(pcchan_gamma,:);
    g_noise2  =  noise_after_gamma(pcchan_gamma,:);
    % snr = signal/denoise
    g_snr1    = g_signal1 ./g_noise1;
    g_snr2    = g_signal2./g_noise2;
    
       
%     %% BROADBAND
%     signal_before_bb = nanmedian(before.w_pwr(~denoisedData.bad_channels,:,:),3);
%     noise_before_bb = nanstd(before.w_pwr(~denoisedData.bad_channels,:,:),[],3);
%     
%     signal_after_bb = nanmedian(after.w_pwr,3);
%     noise_after_bb = nanstd(after.w_pwr,[],3);
%     
%     snr_before_bb = signal_before_bb ./ noise_before_bb;
%     snr_after_bb = signal_after_bb ./ noise_after_bb;
%     
%     % max across 3 conditions 
%     finalsnr = cat(2, snr_before_bb, snr_after_bb);
%     % max across before and after
%     finalsnr = max(finalsnr');
%     % sort
%     [~,idx] = sort(finalsnr,'descend');
%     % find the top 10
%     pcchan_bb = false(size(denoisedData.results.noisepool));
%     pcchan_bb(idx(1:10))= 1;
%     
%     
%     % signal and noise before denoising
%     g_signal1 = signal_before_bb(pcchan_bb,:);
%     g_noise1  = noise_before_bb(pcchan_bb,:);
%     % signal and noise after denoising 
%     g_signal2 = signal_after_bb(pcchan_bb,:);
%     g_noise2  =  noise_after_bb(pcchan_bb,:);
%     % snr = signal/denoise
%     g_snr1    = g_signal1 ./g_noise1;
%     g_snr2    = g_signal2./g_noise2;
    
    condColors = jet(9);
    
    % plot each condition as a different color 
    % signal 
    subject_num2 = subject_num-3;
    subplot(3,3,subject_num2); cla; hold on;
    for nn = 1:9
        plot(g_signal1(nn,:),g_signal2(nn,:),'o','color',condColors(nn,:), 'MarkerSize',8, 'LineWidth',2);
    end
    axis square;
    axismax = max([g_signal1(:);g_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    title(sprintf('S%d : signal', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    
    % noise
    subplot(3,3,subject_num2+length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(g_noise1(nn,:),g_noise2(nn,:),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([g_noise1(:); g_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k','LineWidth',2);
    axis square;
    title(sprintf('S%d : noise', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    % snr
    subplot(3,3,subject_num2+2*length(which_data_to_visualize)); cla; hold on;
    for nn = 1:9
        plot(g_snr1(nn,:),g_snr2(nn,:),'o','color',condColors(nn,:),'MarkerSize',8, 'LineWidth',2);
    end
    axismax = max([g_snr1(:); g_snr2(:)])*1.2;
    xlim([0,axismax]);  set(gca,'FontSize',18); ylim([0,axismax]);  set(gca,'FontSize',18); line([0,axismax],[0,axismax],'color','k', 'LineWidth',2);
    set(gca,'FontSize',18)
    axis square;
    title(sprintf('S%d : SNR', which_data_to_visualize(subject_num2)), 'FontSize',18);
    makeprettyaxes(gca,18,18);
    
    drawnow;
end