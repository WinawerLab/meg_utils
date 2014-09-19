function ssmeg_make_spectra_figures(sensor_data, data_pth, ab, bbPoly, frequencies, amps_on_full, amps_off_full, npcs, denoise_via_pca)

% This function produces a power spectrum plot (for frequencies ranging
% from 10-120 Hz), for every single channel. The signal of both conditions
% are fitted to a line, in order to calculate the broadband signal.

% INPUTS:
% sensor_data       = Amplitude data to plot on head (Difference broad
%                     band/sum broadband).
% data_pth          = Data path in order to safe the figures
% ab                = i.e. Asynchronous Broadband, from ssmeg_broadband
%                     function, which is the broadband based on linear fit
% bbPoly            = The broadband polynomial fits for each electrode
% frequencies       = The frequency range used to fit the data, received from 
%                     ssmeg_broadband function. Only used to plot the data.
% amps_*_full       = Amplitude data for the two conditions separated, for
%                     plotting purposes.
% npcs              = Nr of principle components used in the denoising
%                     analysis. Used for the title of the plots.

% This function needs the fieldtrip toolbox

%% Load data header for positions of the electrodes

% Safe figures in the end?
save_figures = true;

% Get header info
if notDefined('data_hdr')
   data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
end

% If there are no pcs defined, we assume no pcs are used (i.e. nr pcs = 0)
if ~exist('npcs', 'var')
    npcs = 0;
end

% Make new cfg framework for fieldtrip
cfg=[];

% Configure our Field trip layout
cfg.layout          = ft_prepare_layout(cfg, data_hdr);
cfg.update          = 'no';
cfg.style           = 'straight';
cfg.electrodes      = 'on';
cfg.showlabels      = 'markers';
cfg.colorbar        = 'yes';
% cfg.maplimits     = 'maxmin';
cfg.data            = sensor_data';
cfg.emarkersize     = 6;

% Make sure the bad electrodes (with NaN's) are replaced by a value 
%(median in this case)
for ii = 1:length(cfg.data)
    if isnan(cfg.data(ii))
        cfg.data(ii) = nanmedian(sensor_data);
    end
end

%% Loop over data points
for chan = 1:157
    
    % Fill in the electrode of interest
    cfg.highlight = chan;
    
    % Plot data on head
    figure(101); clf; set(gcf, 'Color', 'w', 'Position', [0 0 1000 1000]);
    handle1 = subplot(211);
    position = get(handle1, 'Position');
    set(handle1, 'Position', [position(1) position(2)+.125 position(3)-.1 position(4)-.1]);
    topoplot(cfg,cfg.data(1:157));
    
    
    % Plot poly-line fit with off/on spectra 
    frequencies.all_i = 1:length(frequencies.all);
    handle2 = subplot(212);
    position2 = get(handle2, 'Position');
    set(handle2, 'Position', [position2(1) position2(2)-.05 position2(3) position2(4)+.25]);
    
    % Calculate fitted line for higher frequencies, for every channel, for
    % the two conditions
    mnbb    = polyval(bbPoly(:,chan), log(frequencies.all));
    xon     = log(ab(1,chan)) - polyval(bbPoly(:,chan), log(12));
    xoff    = log(ab(2,chan)) - polyval(bbPoly(:,chan), log(12));
    mnon    = mnbb + xon;
    mnoff   = mnbb + xoff;
    
    % Plot lines for on and off
    plot(frequencies.all, exp(mnon), 'k--', 'LineWidth', 2); hold on;
    plot(frequencies.all, exp(mnoff), '--', 'Color', [.7 .7 .7], 'LineWidth', 2);
    
    % Plot the data
    plot(frequencies.all, nanmedian(amps_on_full(frequencies.all_i,:,chan).^2,2), 'k', 'LineWidth', 2); hold on;
    plot(frequencies.all, nanmedian(amps_off_full(frequencies.all_i,:,chan).^2,2), '-', 'Color', [.7 .7 .7], 'LineWidth', 2); hold on;
 
    % Plot the data points used for the AB linear fit (circles)
    plot(frequencies.ab, nanmedian(amps_on_full(frequencies.ab_i,:,chan).^2,2), 'o', 'Color', 'k', 'MarkerSize', 8);
    plot(frequencies.ab, nanmedian(amps_off_full(frequencies.ab_i,:,chan).^2,2), 'o', 'Color', [.8 .8 .8], 'MarkerSize', 8);
    

    %% Layout
    
    xlabel('Frequency [Hz]', 'FontSize', 16);
    ylabel('Power [pico Tesla]', 'FontSize', 16);
    legend('Fitted line flickering periods', 'Fitted line blank periods', 'Flickering periods', 'Blank periods', 'Location', 'SouthWest' );
    title(sprintf('Power spectrum of channel nr %d, principle component nr %d', chan, npcs), 'FontSize', 18);
    set(gca, 'YScale', 'log', 'XScale', 'log', 'XLim', [10, frequencies.all(end)+5], ...
        'XTick', (1:10) * 12, 'XGrid', 'on', 'Fontsize', 16);
   
    % Save spectrum figure
    if save_figures
        if denoise_via_pca
            filename= ['pc_' int2str(npcs) '_denoised_powerspectrum_channel_' int2str(chan) '.eps'];
            hgexport(101,fullfile(data_pth,'figures/headplots_spectra_denoised',filename));
        else
            filename= ['powerspectrum_channel_' int2str(chan) '.eps'];
            hgexport(101,fullfile(data_pth,'figures/headplots_spectra',filename));
        end
    end
    %% Go on to the next channel
%      waitforbuttonpress;
    

end
    
end




