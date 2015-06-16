function cond_diff_headplot(data, which_conditions, which_model, type)
%% Function to reproduce Figure 9 (Spatialmap) for right minus left denoised activations 
%
% cond_diff_headplot(bb, which_model, type)
%
% INPUTS:
%   data: data structure from denoising function, must include
%       results field, with subfields origmodel and finalmodel, and further
%       subfields beta_md and beta_se.
%   which_conditions: which conditions to compare, must be: right_left,
%       full_right, or full_left)
%   which_model: origmodel or finalmodel
%   type: do you want to compare signal, noise, or SNR: must be 'S', 'N',
%   or 'SNR'
%   
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the broadband signals in the contrast right minus left
% after using the denoising algorithm. 

%%
% Condition numbers refer to:
%   1: Full field (bilateral)
%   2: Left (from stimulus numbers 5 and 6)
%   3: Right (from stimulus numbers 7 and 8) 
%
% See s_SSEEG_Analysis and stimulus files such as onOffLeftRight_600x600params*

switch which_conditions
    case 'left_right'
        first_cond = 2; second_cond = 3;
        comparison = 'Left Minus Right';
    case 'full_right'
        first_cond = 1; second_cond = 3;
        comparison = 'Full Minus Right';
    case 'full_left'
        first_cond = 1; second_cond = 2;
        comparison = 'Full Minus Left';
end 

%% Plot spatial map figures: right minus left stimulation for broadband component after denoising.

    results = data.results;
    
    ab_snr1 = getsignalnoise(results.(which_model), first_cond, type); % Left
    ab_snr2 = getsignalnoise(results.(which_model), second_cond, type); % Right
    
    ab_snr_diff = ab_snr1 - ab_snr2;
    fH = figure;
    plotOnEgi(ab_snr_diff(~data.badChannels)); colorbar;
    set(colorbar, 'Limits', [-2.5 2.5])
    caxis([-2.5 2.5])
    title(sprintf('%s, %d Noise Channels', comparison, sum(data.results.noisepool)));
    
end
