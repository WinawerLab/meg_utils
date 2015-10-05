% s_GAMMA_investigate_spectra
%
% Script to plot (mean) spectra of different images and channels of the
% Gamma experiment. For now, script can only plot data from the 1 second
% image presentation experiment, before denoising (so 1000 time points).




%% Get paths
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';

subject_num                   = 9;
save_images                   = true;

% meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})

d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

condition_names               = {   ...
    'White Noise' ...
    'Binarized White Noise' ...
    'Binarized Pink Noise' ...
    'Binarized Brown Noise' ...
    'Gratings (0.36 cpd)' ...
    'Gratings (0.73 cpd)' ...
    'Gratings (1.45 cpd)' ...
    'Gratings (2.90 cpd)' ...
    'Plaid'...
    'Blank'};

rgb_pink  = [228,43,145]/255;
rgb_green   = [141,200,108]/255;
rgb_purple    = [93,12,139]/255;
rgb_grey    = [75, 75, 75]/255;

color_scheme = [rgb_pink;rgb_pink;rgb_pink;rgb_pink; ...
    rgb_green;rgb_green;rgb_green;rgb_green;rgb_purple];

satValues       = 1-linspace(0.1,1,4);
colorRGB_pink   = varysat(rgb_pink,satValues);
colorRGB_noise  = varysat(rgb_green,satValues);

%% Load data
load_pth    = fullfile(project_pth, subj_pths{subject_num}, 'processed');
datasets    =  dir(fullfile(load_pth, '*boot*'));
before      = load(fullfile(load_pth, datasets(1).name));
after       = load(fullfile(load_pth, datasets(2).name));


denoisedData = dir(fullfile(load_pth, sprintf('s0%d_denoisedData*',subject_num+1)));
denoisedData = load(fullfile(load_pth,denoisedData(1).name));

spectral_data_files  = dir(fullfile(load_pth, 'spectral_data*.mat'));
spectral_data_before = load(fullfile(load_pth,spectral_data_files(1).name));
% spectral_data_after  = load(fullfile(load_pth,spectral_data_files(2).name));

model_fit_before = nanmedian(before.fit_f2,4);
data_before     = nanmedian(spectral_data_before.spectral_data_boots,4);



t = (1:1000)/1000;
f = (0:length(t)-1)/max(t);
for chan =  1:20
        %% Plot mean across conditions
        figure(3), clf;
     
        subplot(211);  hold all;
        % mean noise
        plot(f,smooth(mean(data_before(:,1:4,chan),2),2)', 'LineWidth',2);
%         plot(f,10.^model_fit_before(ii,:,chan),'color', color_scheme(ii,:,:), 'LineWidth',4); hold on;
        
        % mean gratings
        plot(f,smooth(mean(data_before(:,5:8,chan),2),2)', 'LineWidth',2);
%         plot(f,10.^model_fit_before(10,:,chan),'color',rgb_grey,'LineWidth',4);

        % plaid
        plot(f,smooth(data_before(:,9,chan),2)',  'LineWidth',2);

        
        % baseline
        plot(f,smooth(data_before(:,10,chan),2)', 'color',rgb_grey, 'LineWidth',2);
        
        
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 65])
         ylim([5 40])
        set(gca,'XTick',[30:10:80],'XGrid','on')
        box(gca,'off');        set(gcf, 'color','w')
        title('Mean across conditions, within subject and channel', 'FontSize',18)
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend('Noise', 'Gratings','Plaid','Baseline');
%         hgexport(gcf, fullfile(project_pth, subj_pths{subject_num}, 'figs',sprintf('data_modelfit_before_chan%d_meanconds_2',chan)));

        %% Plot selected channels of a subject and look whether gamma oscillation peak differs per condition
        subplot(212); cla; hold all;
        for ii = 1:9;
        plot(f,smooth(mean(data_before(:,ii,chan),2),2)',  'color', color_scheme(ii,:,:), 'LineWidth',2);
%         plot(f,10.^model_fit_before(ii,:,chan),'color', color_scheme(ii,:,:), 'LineWidth',4); hold on;
        end

        
        % baseline
        plot(f,smooth(data_before(:,10,chan),2)', 'color',rgb_grey, 'LineWidth',2);
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([10 65])
        ylim([3 50])
        set(gca,'XTick',[30:10:80],'XGrid','on')
        box(gca,'off');        set(gcf, 'color','w')
        title(num2str(chan), 'FontSize',18)
        
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
%         legend('Noise', 'Gratings ','Plaid','Baseline');

     waitforbuttonpress;
    end
    
model_fit_after = nanmedian(after.fit_f2,4);
data_after     = nanmedian(spectral_data_after.spectral_data_boots,4);

% In case we want to dismiss some frequencies for plotting
%     data_after(29,:,:) = NaN;
%     data_after(32,:,:) = NaN;
%     data_after(34,:,:) = NaN;
%     data_after(49,:,:) = NaN;
%     data_after(97,:,:) = NaN;
%     data_after(102,:,:) = NaN;
%     data_after(141,:,:) = NaN;
%     data_after(150,:,:) = NaN;

t = (1:800)/1000;
f = (0:length(t)-1)/max(t);
for chan = 40:80
    for ii = 1:9
        figure(ii), clf;
        
        
        plot(f,10.^model_fit_after(ii,:,chan),'color', color_scheme(ii,:,:), 'LineWidth',4); hold on;
        plot(f,data_after(:,ii,chan), 'color', color_scheme(ii,:,:), 'LineWidth',2);
        
        plot(f,10.^model_fit_after(10,:,chan),'color',rgb_grey,'LineWidth',4);
        plot(f,data_after(:,10,chan), 'color', rgb_grey, 'LineWidth',2);
        
        %                         hgexport(gcf, fullfile(project_pth, subj_pths{subject_num}, 'figs',sprintf('data_modelfit_after_chan%d_allcond',chan,ii)));
        
        set(gca, 'YScale','log','XScale','log','LineWidth',2)
        xlim([30 200])
        ylim([10.^.65, 10.^1.2])
        box(gca,'off');        set(gcf, 'color','w')
        %             title(condition_names{ii}, 'FontSize',18)
        xlabel('Frequency (Hz)','FontSize',18)
        ylabel('Power','FontSize',18)
        legend(condition_names{ii},'Data','Baseline');
        
        hgexport(gcf, fullfile(project_pth, subj_pths{subject_num}, 'figs',sprintf('data_modelfit_after_chan%d_cond%d_2',chan,ii)));
        
        
    end
    
    
    
end









