% s_GAMMA_conditions_spectra
% 
% Plots the spectra of each condition for every MEG channel

%% Paths and options
project_pth                    = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                       = '*_Gamma_*subj*';

% options
subjects                       = 20;
save_images                    = true;

dateFormat   = 'mm_dd_yy';
suffix       = 'localregression_multi_100_boots'; 
figSuffix    = sprintf('%s_%s' , datestr(now, dateFormat), suffix);

d                = dir(fullfile(project_pth, data_pth));
subj_pths        = struct2cell(d);
subj_pths        = subj_pths(1,:);

colormap         = parula(10);
scrsz            = get(0,'ScreenSize');

%% Loop over subjects
for this_subject = subjects
   condition_names = gamma_get_condition_names(this_subject);
   
   % this_subject -1
   load_pth = fullfile(project_pth, subj_pths{this_subject -1}, 'processed');
   thisFile = sprintf('spectral_data_%s.mat', suffix);
   spectral_data_files  = dir(fullfile(load_pth, thisFile));
   spectral_data        = load(fullfile(load_pth, spectral_data_files(1).name));
   
   data = spectral_data.spectral_data_boots;
   
   if numel(size(data)) == 4
       data = nanmedian(spectral_data.spectral_data_boots, 4);
   end
   
   num_timepoints = size(data, 1);
   fs = 1000;
   t = (1:num_timepoints)/fs;
   f = (0:length(t)-1)/max(t);
   
   % change the frequency range to which the data is fit, excludes noise
   % frequency bands
   f_use4fit = f((f< 57) | (f>=62 & f <= 117) | (f>=123 & f <= 178) | (f>=182 & f <= 200));
   f_sel = ismember(f,f_use4fit); % or use f_sel from the data files
   
   %% loop across channels
   for chan = 1:157
      % plot all conditions on the same spectrogram
      figure(1), clf;
      set(gcf, 'color', 'w');
      hold all;
      
      for ii = 1:length(condition_names)
          plot(f(f_sel), smooth(data(f_sel, ii, chan),2)', 'color', colormap(ii,:,:), 'LineWidth', 2);
      
      end
   end
   
end