function makeRLheadplot()
%% Function to reproduce Figure 9 (Spatialmap) for right minus left denoised activations 
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the broadband signals in the contrast right minus left
% after using the denoising algorithm. 
%

%% STILL UNDER CONSTRUCTION

%% Plot spatial map figures: right minus left stimulation for broadband component after denoising.
figure('position',[1,600,1400,800]);
whichmodel = 'finalmodel';
for k = 1:length(whichSubjects)

    data = dataAll{k};
    results = data{1}.results;
    
    subplot(2,4,k);  
    ab_snr1 = getsignalnoise(results.(whichmodel), 2, 'SNR'); % Left
    ab_snr2 = getsignalnoise(results.(whichmodel), 3, 'SNR'); % Right
    ab_snr_diff = to157chan(ab_snr2-ab_snr1,~data{1}.badChannels,'nans');
    
    if dataType == MEG    
    [~,ch] = megPlotMap(ab_snr_diff,[-5,5],gcf,jmaColors('coolhotcortex'));
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-5:1:5);
    makeprettyaxes(ch,9,9);
    
    elseif dataType == EEG
    plotOnEgi(ab_snr_diff(~data{1}.badChannels));
    
end        
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure9ab_bbRightMLeft_after'),[],0,'.',1);
end
