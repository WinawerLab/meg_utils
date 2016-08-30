% Visualize spectra and gamma + bb model fits


% Note about the different datasets one can choose from:
%                   Peak restr. Freq used
% Modelfit 1:       35-80 Hz    (f>=35 & f < 40)  | (f > 40 & f <= 57) | (f>=65 & f <= 115)  | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 2:       40-80 Hz    (f>=35 & f <= 58) | (f>=62 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 3:       40-80 Hz    (f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 4:       42-80 Hz    (f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 5:       46-80 Hz    (f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));

% Choose project path
project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
% project_pth                   = '~/matlab/git/meg_utils/experiments/gamma_meg/HPC/Data/';
data_pth                      = '*_Gamma_*subj*';

% Define parameters
fs                            = 1000;
which_session_to_visualize    = 15; %8;%[6,7,8,9,10,11,12,14,15,16]
save_images                   = true;

% Get the folders of the subjects in the Gamma experiment
% Note: the session numbers are always subject_nr +1
d = dir(fullfile(project_pth, data_pth));
subj_pths = struct2cell(d);
subj_pths = subj_pths(1,:);

% Define colors
rgb_pink  = [228,43,145]/255;
rgb_green   = [141,200,108]/255;
rgb_purple    = [93,12,139]/255;
rgb_grey    = [75, 75, 75]/255;

color_scheme = [rgb_pink;rgb_pink;rgb_pink;rgb_pink; ...
    rgb_green;rgb_green;rgb_green;rgb_green;rgb_purple];

% Define saturated colors in case you want to plot all fits in one figure
satValues       = 1-linspace(0.1,1,4);

% colorRGB_pink   = varysat(rgb_pink,satValues);
% colorRGB_noise  = varysat(rgb_green,satValues);

%% Loop over data sets
for session_num = which_session_to_visualize
    
    % Get condition names (different for the first and last 6 subjects)
    [condition_names, baseline_condition]  = gamma_get_condition_names(session_num);
    
    fprintf('Session number %d of %d\n', ...
        find(which_session_to_visualize==session_num), ...
        length(which_session_to_visualize)); drawnow();
    
    %% Define paths and load data
    load_pth    = fullfile(meg_gamma_get_path(session_num), 'processed');
    
    % Go to specific subject folder and find datasets
    datasets    = dir(fullfile(load_pth, '*local*'));
    
    % Define which datasets you want to load
    for ii = 1:numel(datasets)
        fprintf('%d: %s\n', ii, datasets(ii).name)
    end
    
    % Load the parameters of the model fit (Gaussian weight, broadband
    % weight, etc)
    %   Before is without denoising, after is with denoising
    before_dataset = input('Modelfit number before denoising?');
    after_dataset  = input('Modelfit number after denoising?');
    
    before  = load(fullfile(load_pth, datasets(before_dataset).name));
    after   = load(fullfile(load_pth, datasets(after_dataset).name));
    
    % Get denoised dataset in case we need to define the badChannels
    denoisedData = dir(fullfile(load_pth, sprintf('s%02d_denoisedData*',session_num)));
    denoisedData = load(fullfile(load_pth,denoisedData(1).name));
    
    % Get spectra for before and after denoising
    spectral_data_files  = dir(fullfile(load_pth, '*spectral_data_*local*.mat'));
    
    for ii = 1:numel(spectral_data_files)
        fprintf('%d: %s\n', ii, spectral_data_files(ii).name)
    end
    before_dataset = input('Get spectra without denoising from which results file number?');
    after_dataset  = input('Get spectra parameters with denoising from which results file number?');
    
    spectral_data_before = load(fullfile(load_pth,spectral_data_files(before_dataset).name));
    spectral_data_after  = load(fullfile(load_pth,spectral_data_files(after_dataset).name));
    
    % Take the median across bootstraps for modelfit and spectra
    model_fit_before = nanmedian(before.fit_f2,4);
    %     model_fit_before = before.fit_f2;
    data_before      = nanmedian(spectral_data_before.spectral_data_boots,4);
    
    % Do the same but then for model and spectra after denoising
    model_fit_after = nanmedian(after.fit_f2.^2,4);
    data_after      = nanmedian(spectral_data_after.spectral_data_boots,4);
    
    
    % Define time and frequencies to use
     
    num_channels = size(spectral_data_before.spectral_data_boots,3);
    
    
    %% Plot Spectra
    
    %     % plots all conditions and model fit of every channel
    %     figure(1);clf
    %     for chan = 1:157;
    %         plot(...
    %             f, exp(mean(log(spectral_data_before.spectral_data_boots(:,:,chan,:)),4)), 'k', ...
    %             f, exp(model_fit_before(:,:,chan))', 'r')
    %         set(gca, 'YScale', 'log','XScale', 'log', 'XLim', [10 200])
    %         title(chan)
    %         waitforbuttonpress;
    %     end
    
    
    %% Plot
    fH = figure('position', [1,400,400,800]);
    xt = [35 50 100 150 200];
    
    for before_or_after = 1:2
        switch before_or_after
            case 1
                data = data_before;
                model_fit = model_fit_before;
                params = before;
                str = 'before';
                t = (1:1000)/1000;
                f = (0:length(t)-1)/max(t);
                f_sel = intersect(f, before.f_use4fit);
                
            case 2
                data = data_after;
                model_fit = model_fit_after;
                str = 'after';
                params = after;
                t = (1:800)/1000;
                f = (0:length(t)-1)/max(t);
                f_sel     = intersect(f, after.f_use4fit);
                    
        end
        
        for chan = 23;%:num_channels
            set(fH, 'name', sprintf('Channel %d', chan));
            for ii = 1:9
                
                clf;
                %             subplot(10,1,ii)
                
                %             plot(f(f_sel),10.^model_fit_before(ii,f_sel,chan), '-o','color', color_scheme(ii,:,:), 'LineWidth',2); hold on;
                
                plot(f,exp(model_fit(ii,:,chan)), '-','color', color_scheme(ii,:,:), 'LineWidth',2); hold on;
                
                plot(f,data(:,ii,chan), 'color', color_scheme(ii,:,:), 'LineWidth',1);
                %             plot(f,mean(data_before(:,ii,chan),2), 'color', color_scheme(ii,:,:), 'LineWidth',2);
                
                
                %             fit_withoutgamma = model_fit_before(ii,f_sel,chan) - ...
                %                 before.w_gauss_mn(chan, ii) * ...
                %                 0.04*sqrt(2*pi)*normpdf(log(f_use4fit),before.gauss_f(chan, ii),0.04);
                
                plot(f,exp(model_fit(baseline_condition,:,chan)),'color',rgb_grey,'LineWidth',2);
                
                %             plot(f(f_sel),10.^model_fit_before(10,f_sel,chan),'color',rgb_grey,'LineWidth',4);
                plot(f,data(:,baseline_condition,chan), 'color', rgb_grey, 'LineWidth',1);
                
                
                set(gca, 'YScale','log','XScale','log','LineWidth',1)
                xlim([35 200])
                set(gca,'XTick',[25 35 50 75 100 150 200])
                set(gca,'XTickLabel',{'25','35','50','75','100','150','200'})
%                 set(gca,'XLim', [35 200]) %200axis([xt(1) yl])
%                 set(gca,'YLim', [10.^-1 10.^3]) %200axis([xt(1) yl])
                set(gca,'YLim', [10.^0 10.^2.2]) % The right axes for MEG plot alone
%                 set(gca,'YLim',[10.^-2 10.^2.25]) % The axis of ECOG

                
                set(gca,'XGrid','on')
%                 set(gca,'YTick',10^-2:50:10^2.25)
                set(gca,'box', 'off');        set(gcf, 'color','w')
                makeprettyaxes(gca,9,9)
                
                title(...
                    sprintf('%s, Broadband: %5.3f, Gamma: %5.3f at %3.1f Hz', ...
                    condition_names{ii}, params.w_pwr_mn(chan, ii)-params.w_pwr_mn(chan, baseline_condition), ...
                    params.w_gauss_mn(chan, ii)-params.w_gauss_mn(chan, baseline_condition),...
                    exp(params.gauss_f(chan, ii))),...
                    'FontSize', 20)
                
                xlabel('Frequency (Hz)','FontSize',20,'FontWeight','bold')
                ylabel('Power (pico Tesla)','FontSize',20,'FontWeight','bold')
                legend({sprintf('Modelfit %s', condition_names{ii}), sprintf('Data %s', condition_names{ii}), 'Modelfit Baseline', 'Data Baseline'},'FontSize',20,'FontWeight','bold');
                legend('boxoff')
                if save_images
                    poster_path = sprintf('/Users/winawerlab/Google Drive/FYP/15_Gamma_9_29_2015_subj028/MEG_AXES', session_num);
                    
                    hgexport(gcf, fullfile(poster_path,...
                        sprintf('data_modelfit_%s_chan%d_cond%d_local_regression_boot100', str, chan,ii)));
                end
                %             waitforbuttonpress
            end
            
        end
    end
    
  

end