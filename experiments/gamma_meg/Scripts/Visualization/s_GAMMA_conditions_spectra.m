% s_GAMMA_conditions_spectra
% 
% Plots the spectra of each condition for every MEG channel

%% Paths and options
project_pth                    = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                       = '*_Gamma_*subj*';

% options
subjects                       = 23;
save_images                    = true;

dateFormat   = 'mm_dd_yy';
suffix       = '100boots'; 
figSuffix    = sprintf('%s_%s' , datestr(now, dateFormat), suffix);

d                = dir(fullfile(project_pth, data_pth));
subj_pths        = struct2cell(d);
subj_pths        = subj_pths(1,:);

colormap         = parula(13);
scrsz            = get(0,'ScreenSize');

%% Loop over subjects
for this_subject = subjects
   condition_names = gamma_get_condition_names(this_subject);
   
   % this_subject -1
   load_pth = fullfile(project_pth, subj_pths{this_subject -1}, 'processed');
   thisFile = sprintf('*spectralData_%s*.mat', suffix);
   spectral_data_files  = dir(fullfile(load_pth, thisFile));
   spectral_data        = load(fullfile(load_pth, spectral_data_files(1).name));
   
   data = spectral_data.spectralDataBoots;
   
   if numel(size(data)) == 4
       data = nanmedian(spectral_data.spectralDataBoots, 4);
   end
   
   num_timepoints = size(data, 1);
   fs = 1000;
   t = (1:num_timepoints)/fs;
   f = (0:length(t)-1)/max(t);
   
   % change the frequency range to which the data is fit, excludes noise
   % frequency bands
%    f_use4fit = f((f< 57) | (f>=62 & f <= 117) | (f>=123 & f <= 178) | (f>=182 & f <= 200));
   f_use4fit = f;
   f_sel = ismember(f,f_use4fit); % or use f_sel from the data files
   
   %% loop across channels
   for chan = 1:size(data,3);
      % plot all conditions on the same spectrogram
      figure(1), clf;
      set(gcf, 'color', 'w');
      
      
      for ii = 1:length(condition_names)
          
          subplot(151);
          hold all;
%           plot(f(f_sel), smooth(data(f_sel, ii, chan),2)', 'color', colormap(ii,:,:), 'LineWidth', 2);
          plot(f, data(f_sel, ii, chan)', 'color', colormap(ii,:,:), 'LineWidth', 2);
          set(gca,'XScale', 'log', 'YScale', 'log', 'Xlim',[1 200]);
          xlabel('Frequency (Hz)'); ylabel('Power (fT)');
          legend(condition_names{:}); legend('Location', 'best');
          
          subplot(152);
          hold all;
          plot(f, data(f_sel, ii, chan)', 'color', colormap(ii,:,:), 'LineWidth', 2);
          plot([f(30), f(30)]+1, [10.^-29, 10.^-25], 'k--')
          plot([f(50), f(50)]+1, [10.^-29, 10.^-25], 'k--')
          plot([f(58), f(58)]+1, [10.^-29, 10.^-25], 'k--')

          set(gca,'XScale', 'log', 'YScale', 'log','Xlim',[20 80])
          xlabel('Frequency (Hz)'); ylabel('Power (fT)');
          
          

      end
      
      
      subplot(153)
      megPlotMap(squeeze(mean(data(f(30)+1, :, :),2)),[0 max(squeeze(mean(data(f(30)+1, :, :),2)))],[],[],'30 Hz topography')
      subplot(154)
      megPlotMap(squeeze(mean(data(f(50)+1, :, :),2)),[0 max(squeeze(mean(data(f(50)+1, :, :),2)))],[],[],'50 Hz topography')
      subplot(155)
      megPlotMap(squeeze(mean(data(f(58)+1, :, :),2)),[0 max(squeeze(mean(data(f(58)+1, :, :),2)))],[],[],'58 Hz topography')
      
      hgexport(1,sprintf('~/Desktop/spectra_chan%d',chan));
     
      
   end
   
  
   
   
end