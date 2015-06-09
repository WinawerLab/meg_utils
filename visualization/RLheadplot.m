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

switch which_conditions
    case right_left
        first_cond = 2; second_cond = 3;
            comparison = 'Right Minus Left';
    case full_right
        first_cond = 1; second_cond = 2;
            comparison = 'Full Minus Right';
    case full_left
        first_cond = 1; second_cond = 3;
            comparison = 'Full Minus Left';
        
%% Plot spatial map figures: right minus left stimulation for broadband component after denoising.

    results = data.results;
    
    ab_snr1 = getsignalnoise(results.(which_model), first_cond, type); % Right
    ab_snr2 = getsignalnoise(results.(which_model), second_cond, type); % Left
    
    ab_snr_diff = ab_snr1 - ab_snr2;
    plotOnEgi(ab_snr_diff(~data.badChannels));
    
    title(sprintf(first_cond))
    
end
