function ts_denoised = meg_environmental_denoising(ts,...
    environmental_channels,data_channels, produce_figures, save_data)

%% Description of function

% This function uses the timeseries of the three MEG noise channels
% (sitting on top of the other 157, so they get environmental/external
% noise). The timeseries will be used in a linear regression, where the
% residuals will be used as the new 'denoised' data. It will make some
% before and after denoising plots.
%
% INPUTS:
%  ts                        MEG timeseries (number of timepoints by number 
%                                of epochs  by number of channels)
%  environmental_channels    vector of channel numbers that are
%                               nonphysiological. These will be used to
%                               regress out noise from all other channels.
%  data_channels             vector of channel numbers to denoise
%  produce_figures           boolean. If true, make some plots to compare
%                               pre and post denoising
%
% OUTPUTS:
%  ts                        denoised time series

%% Deal with inputs
if ~exist('produce_figures', 'var') || isempty(produce_figures), produce_figures = 0; end
if ~exist('save_data', 'var')       || isempty(save_data),       save_data = 0;       end

%% Define timeseries of conditions


%% Make empty arrays for regressed 'clean' data
ts_denoised = ts;

% Start regression, keep residuals
warning off stats:regress:RankDefDesignMat
for channel = data_channels; 
    fprintf('[%s]: Channel %d\n', mfilename, channel); 
    for epoch = 1:size(ts,2); % Epoch size is the same for every condition (i.e. 180 except for session 3 (=168))
        
        %%% ON PERIODS %%%
        
        % Full
        [~,~,R] = regress(ts(:,epoch,channel),[squeeze(ts(:,epoch,environmental_channels)) ones(size(ts,1),1) ]);
        ts_denoised(:,epoch,channel) = R;
        
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

%% For debugging: Make figures of all the raw epochs 
if produce_figures

    % And for two visual channels
    chan_1 = 1;
    figure; plot(squeeze(mean(ts_denoised(:,:,chan_1),2)),'r'); hold on;
    plot(squeeze(mean(ts(:,:,chan_1),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (Picotesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', chan_1))
    legend('Denoised','Raw')


    chan_1 = 14;
    figure; plot(squeeze(mean(ts_denoised(:,:,chan_1),2)),'r'); hold on;
    plot(squeeze(mean(ts(:,:,chan_1),2)),'b')
    xlabel('Time (ms)')
    ylabel('Amplitude (Picotesla)')
    title(sprintf('Before and after denoising - Timeseries of channel nr %d', chan_1))
    legend('Denoised','Raw')


end


