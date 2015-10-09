% s_GAMMA_test_fits
% loads a 1 channel timeseries to test different fitting situations


%% parameters

project_pth                     = '/Volumes/server/Projects/MEG/Gamma/Data';

fs                            = 1000; 
epoch_start_end               = [0.050 1.049];
blank_condition               = 10;
intertrial_trigger_num        = 11;    
subject                       = 6;
channel                       = 13;


meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'});

condition_names  = gamma_get_condition_names(subject);
num_conditions   = length(condition_names);


%% load ts 

load(fullfile(project_pth,sprintf('ts_subj%d_ch%d', subject, channel)));

ts = test_ts;

conditions_unique = unique(conditions);

num_epochs = size(ts,2);

t = (1:size(ts,1)); %timepoints


%% -------Spectral Analysis--------
% case 1: fourier transform without smoothing 
% case 2: fourier transform using sliding window
% case 3: fourier transform using sliding hanning window

%% case 1
% t = (1:size(ts,1))/fs;
% f = (0:length(t)-1)/max(t);
% 
% 
% % fourier transform data
% spectral_data = abs(fft(ts))/length(t)*2;

% take mean across conditions
% spectral_data_mean = zeros(size(ts,1), length(conditions_unique));
% 
% for ii = 1:size(conditions_unique,1);
%     
%     these_epochs = conditions == conditions_unique(ii);
%     
%     these_data = spectral_data(:,these_epochs);
%     
%     spectral_data_mean(:,ii) = exp(nanmean(log(these_data),2));
% %     spectral_data_mean(:,ii) = nanmean(log10(these_data),2);
% 
% end

%% case 2 and 3

use_hanning_window = true;

num_windows = 7; %number of moving windows


N = max(t); % 1000ms per epoch
f = 0:N-1; % frequencies
window_size = (2*N)/(num_windows+1); % size of moving window (ms)
h = (hann(window_size)*2);
h = repmat(h,1,num_epochs);

F_all = abs(fft(ts))/N*2; % fft of 1 window/epoch for comparison
y_windowed = [];
F_windowed = zeros(num_windows,window_size,num_epochs); 

% loop over number of windows
for ii = 1:num_windows
    % isolate ith window of each epoch
    y_windowed = ts((1:window_size)+window_size/2*(ii-1),:);
    
    % multiply each window by hanning function
    if use_hanning_window
        y_windowed = y_windowed .* h;
    end
    
    % fft ith window
    F_windowed(ii,:,:) = abs(fft(y_windowed))./(window_size*2);
end

% take the mean across windows/squeeze
F_windowed = squeeze(mean(F_windowed,1));
% adjusted frequency scale
f_windowed = 0:N/window_size:N-1;

% take mean across conditions
F_windowed_mean = zeros(window_size, length(conditions_unique));

for ii = 1:size(conditions_unique,1);
    
    these_epochs = conditions == conditions_unique(ii);
    
    these_data = F_windowed(:,these_epochs);
    
    F_windowed_mean(:,ii) = exp(nanmean(log(these_data),2));
%     spectral_data_mean(:,ii) = nanmean(log10(these_data),2);
    
    F_all(:,ii) = exp(nanmean(log(F_all(:,these_epochs)),2));

end



for iii = 1:num_conditions
    
    title(condition_names(iii));
    figure(iii); plot(f, F_all(:,iii), 'r', f_windowed, F_windowed_mean(:,iii), 'b')
    set(gca, 'Yscale', 'log', 'Xscale', 'log')
    
    
    waitforbuttonpress
end


    
%% Model Fit Parameters 

%                   Peak restr. Freq used
% Modelfit 1:       35-80 Hz    (f>=35 & f < 40)  | (f > 40 & f <= 57) | (f>=65 & f <= 115)  | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 2:       40-80 Hz    (f>=35 & f <= 58) | (f>=62 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 3:       40-80 Hz    (f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
% Modelfit 4:       42-80 Hz    (f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));


f_use4fit = f((f>=35 & f <= 57) | (f>=63 & f <= 115) | (f>=126 & f <= 175) | (f>=186 & f <= 200));
f_sel=ismember(f,f_use4fit);
num_time_points = round((epoch_start_end(2)-epoch_start_end(1)+0.001)*fs);

out_exp = NaN(num_conditions);     % slope of spectrum in log/log space
w_pwr   = NaN(num_conditions);     % broadband power
w_gauss = NaN(num_conditions);     % gaussian height
gauss_f = NaN(num_conditions);     % gaussian peak frequency
fit_f2  = NaN(num_conditions,num_time_points); % fitted spectrum

warning off 'MATLAB:subsassigndimmismatch'

for cond = 1:num_conditions-1
    
    data_to_fit  = spectral_data_mean(:,cond);
    data_for_baseline = spectral_data_mean(:,blank_condition);
    
    try
    [out_exp(cond), ...
        w_pwr(cond), ...
        w_gauss(cond),...
        gauss_f(cond),...
        fit_f2(cond,:)] = ...
        gamma_fit_data(f,f_use4fit,data_for_baseline,data_to_fit);
    catch ME
         warning(ME.identifier, ME.message)
    end



end

%% Plot spectra

% Define colors
rgb_pink  = [228,43,145]/255;
rgb_green   = [141,200,108]/255;
rgb_purple    = [93,12,139]/255;
rgb_grey    = [75, 75, 75]/255;

color_scheme = [rgb_pink;rgb_pink;rgb_pink;rgb_pink; ...
                rgb_green;rgb_green;rgb_green;rgb_green;rgb_purple];



model_fit = fit_f2;


for ii = 1:9
    %figure(ii),
    clf;
    %             subplot(10,1,ii)
    
    plot(f(f_sel),10.^model_fit(ii,f_sel),'color', color_scheme(ii,:,:), 'LineWidth',4); hold on; %model of all conditions
    plot(f,spectral_data_mean(:,ii), 'color', color_scheme(ii,:,:), 'LineWidth',2);
    %plot(f,mean(ts(:,5:8),2), '-o', 'color', color_scheme(ii,:,:), 'LineWidth',2); %ts gratings
    plot(f(f_sel),10.^model_fit(10,f_sel),'color',rgb_grey,'LineWidth',4); %baseline
    plot(f,spectral_data_mean(:,10), 'color', rgb_grey, 'LineWidth',2);%ts baseline
    
    
    set(gca, 'YScale','log','XScale','log','LineWidth',2)
    xlim([20 200])
    %ylim([3 35])
    set(gca,'XTick',[30:10:80],'XGrid','on')
    box(gca,'off');        set(gcf, 'color','w')
    title(condition_names{ii}, 'FontSize',18)
    xlabel('Frequency (Hz)','FontSize',18)
    ylabel('Power','FontSize',18)
    legend(condition_names{ii}, 'Data', 'Baseline');
    %     %hgexport(gcf, fullfile(project_pth, subj_pths{subject_num}, 'figs',sprintf('data_modelfit_before_chan%d_cond%d_4',chan,ii)));
    waitforbuttonpress
end












