% s_GAMMA_MEG_summarize

project_pth                   = '/Volumes/server/Projects/MEG/Gamma/Data';
data_pth                      = '*_Gamma_*subj*';
% which_session_to_visualize    = [10:12,14:16];%5:9; %
which_session_to_visualize    = 5:9; %
save_images                   = true;
using_denoised_data           = true;
% suffix                        = 'localregression_multi_100';
suffix                        = 'denoisedData_bootstrapped_localregression_multi_100';
%% Derived
fs                            = 1000;

if isempty(which('ft_prepare_layout'))
    meg_add_fieldtrip_paths('/Volumes/server/Projects/MEG/code/fieldtrip',{'yokogawa', 'sqdproject'})
end


condition_names = gamma_get_condition_names(which_session_to_visualize(1));

%% Calculating SNR contrasts

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

snr_w_gauss_summary = zeros(num_contrasts, 157, length(which_session_to_visualize));
snr_w_pwr_summary = zeros(num_contrasts, 157, length(which_session_to_visualize));

%% Loop over data sets
for session_num = which_session_to_visualize
    
    fprintf('Session number %d of %d\n', ...
        find(which_session_to_visualize==session_num), ...
        length(which_session_to_visualize)); drawnow();
    
    %% Define paths and load data
    load_pth    = fullfile(meg_gamma_get_path(session_num), 'processed');
    
    if using_denoised_data,  badChannels = [];
    else badChannels = zeros(1,157); end
    
    results = load(fullfile(load_pth, sprintf('s%02d_%s', session_num, suffix)));
    
    if results.nboot > 1,
        summary_stat = @(x) nanmean(x,3) ./ nanstd(x, [], 3);
    else
        summary_stat = @(x) nanmean(x,3);dbquit
    end
    
    num_conditions = size(results.w_pwr,2);
    num_channels   = size(results.w_pwr,1);
    
    if isempty(badChannels)
        denoisedData   = load(fullfile(load_pth,sprintf('s%02d_denoisedData.mat',session_num)));
        badChannels    = denoisedData.bad_channels;
    end
        
    % compute SNR
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
           
    snr_w_gauss_summary(:,:,which_session_to_visualize == session_num) = ...
        to157chan(snr_w_gauss', ~badChannels,0);
    
    snr_w_pwr_summary(:,:,which_session_to_visualize == session_num) = ...
        to157chan(snr_w_pwr', ~badChannels,0);
    
    
end

%% Plot
save_pth = fullfile(project_pth,'Images');

%% SNR Mesh for Gaussian bump

for data_type = 1:2
    scrsz = get(0,'ScreenSize');
    
    switch data_type 
        case 1
            data = mean(snr_w_gauss_summary(1:9,:,:),3);
            str = 'Gamma';
        case 2
            data = mean(snr_w_pwr_summary(1:9,:,:),3);
            str = 'Broadband';
    end
    
    % gaussian weight for each stimuli
    fH = figure; clf; set(fH, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);  
    set(fH, 'name', sprintf('%s SNR', str))
    plot_range = [-1 1] * ceil(max(abs(data(:))));
    threshold =  1/3 * plot_range(2);
    for c = 1:9
        subplot(3,3,c)
        data_to_plot = data(c,:);
        data_to_plot(abs(data_to_plot) < threshold) = 0;
        ft_plotOnMesh(data_to_plot, contrastnames{c});
        set(gca, 'CLim', plot_range);
        colormap(jmaColors('coolhotcortex'))
    end
    
    if save_images
        if ~exist(save_pth, 'dir'), mkdir(save_pth); end
        if using_denoised_data; postFix = 'denoised'; else postFix = []; end;
        hgexport(fH, fullfile(save_pth,sprintf('Group_avg_Per_Condition_%s_SNR_%s_%s.eps',str, suffix, postFix)));
    end
    
end
