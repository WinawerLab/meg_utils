% s_GAMMA_noise_vs_gratings_spectra
%
% Script to plot spectra of grating and noise images for all channgels channels of the
% Gamma experiment.

%% Get paths
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';

subjects                       = [1:11,13:15];
save_images                   = true;

% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

d = dir(fullfile(project_pth, data_pth));
subj_pths        = struct2cell(d);
subj_pths        = subj_pths(1,:);



colormap         = parula(10);
scrsz            = get(0,'ScreenSize');

for subject = subjects
    condition_names  = gamma_get_condition_names(subject);
    
    %% Load data
    load_pth    = fullfile(project_pth, subj_pths{subject}, 'processed');
    
    spectral_data_files  = dir(fullfile(load_pth, 'spectral_data.mat'));
    spectral_data        = load(fullfile(load_pth, spectral_data_files(1).name));
    
    % Load data that is the mean spectral data per condition (freqs x conds x chan => 1000 x 10 x 157)
    data          = spectral_data.spectral_data_boots;
    
    % Or if data is bootstrapped, take the median across bootstraps
    if numel(size(data)) == 4
        data      = nanmedian(spectral_data.spectral_data_boots,4);
    end
   
    
    % Check number of timepoints (should be 500 for first 3 sessions, 1000 for
    % the others)
    num_timepoints = size(data,1);
    fs = 1000;
    
    t = (1:num_timepoints)/fs;
    f = (0:length(t)-1)/max(t);
    
    f_use4fit = f((f< 57) | (f>=62 & f <= 117) | (f>=123 & f <= 178) | (f>=182 & f <= 200));
    f_sel=ismember(f,f_use4fit);
    
    for chan =  1:157
        %% Plot mean across conditions
        
        figure(1), clf; set(gcf, 'color','w', 'position',[scrsz(1),scrsz(2)/2,scrsz(3)/2,scrsz(4)/2]); 
        
        subplot(1,2,1); cla; hold all;
        
        % noise
        for ii = 1:4
            plot(f(f_sel),smooth(data(f_sel,ii,chan),2)','color',colormap(ii,:,:), 'LineWidth',2);
            
        end
        
        
        % image baseline
        plot(f(f_sel),smooth(data(f_sel,10,chan),2)','o-', 'color','k', 'LineWidth',2);
        
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 150])
        ylim([0 50])
        set(gca,'XTick',[10:10:100],'XGrid','on')
        box(gca,'off');
        title(sprintf('Spectra of noise conditions, subject %d and channel %d',subject,chan), 'FontSize',18)
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend(condition_names{1:4},'Baseline')
        
        
        subplot(1,2,2); cla; hold all;
        % gratings
        for ii = 5:9
            plot(f(f_sel),smooth(data(f_sel,ii,chan),2)','color',colormap(ii,:,:), 'LineWidth',2);
         
        end
        % image baseline
        plot(f(f_sel),smooth(data(f_sel,10,chan),2)', 'o-','color','k', 'LineWidth',2);
        
        
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 150])
        ylim([0 50])
        set(gca,'XTick',[10:10:100],'XGrid','on')
        box(gca,'off');
        title(sprintf('Spectra of grating conditions and plaid, subject %d and channel %d',subject,chan), 'FontSize',18)
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend(condition_names{5:9},'Baseline')
        if save_images
            if ~exist(fullfile(project_pth, subj_pths{subject}, 'figs'),'dir'); mkdir(fullfile(project_pth, subj_pths{subject}, 'figs')); end
            hgexport(gcf, fullfile(project_pth, subj_pths{subject}, 'figs',sprintf('spectral_data_chan%d',chan)));
        end
    end
    
end







