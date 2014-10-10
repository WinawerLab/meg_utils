function ts_pre = meg_environmental_denoising(ts_pre, produce_figures, save_data, verbose)

%% Description of function

% This function uses the timeseries of the three MEG noise channels
% (sitting on top of the other 157, so they get environmental/external
% noise). The timeseries will be used in a linear regression, where the
% residuals will be used as the new 'denoised' data. It will make some
% before and after denoising plots.
%
% INPUTS:
% ts_pre  = concatenated timeseries of all epochs (number of timepoints by number of epochs by number of channels)
%
% OUTPUTS:
% ts_pre  = concatenated timeseries of all the denoised epochs

%% Deal with inputs
if nargin < 2 || isempty(produce_figures)
    produce_figures = 0;
end
if nargin < 3 || isempty(save_data)
    save_data = 0;
end
if nargin < 4 || isempty(verbose)
    verbose = 0;
end

%% Define timeseries of conditions


%% Make empty arrays for regressed 'clean' data
ts       = zeros(size(ts_pre));
noise_channels = [158 159 160];

% Start regression, keep residuals
warning off stats:regress:RankDefDesignMat
for channel = 1:157;
    if verbose
        fprintf('[%s]: Channel %d\n', mfilename, channel);
    end
    for epoch = 1:size(ts_pre,2); % Epoch size is the same for every condition (i.e. 180 except for session 3 (=168))
        
        %%% ON PERIODS %%%
        
        % Full
        [~,~,R] = regress(ts_pre(:,epoch,channel),[squeeze(ts_pre(:,epoch,noise_channels)) ones(size(ts_pre,1),1) ]);
        ts(:,epoch,channel) = R;
        
        clear R
        
    end
end
warning on stats:regress:RankDefDesignMat

%% Save denoised data
if save_data
    fprintf('[%s]: Save data matrix to folder..', mfilename);
    
    save denoised_with_nuissance_data.mat ts
    
    fprintf('[%s]: Done! Data matrix saved to folder..', mfilename);
end

%% Make figures of all the raw epochs 
if produce_figures

    % And for two visual channels
    chan_1 = 1;
    figure(104); plot(squeeze(mean(ts(:,:,chan_1),2)),'r'); hold on;
    plot(squeeze(mean(ts_pre(:,:,chan_1),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (Picotesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', chan_1))
    legend('Denoised','Raw')


    chan_1 = 14;
    figure(105); plot(squeeze(mean(ts(:,:,chan_1),2)),'r'); hold on;
    plot(squeeze(mean(ts_pre(:,:,chan_1),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (Picotesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', chan_1))
    legend('Denoised','Raw')


end


