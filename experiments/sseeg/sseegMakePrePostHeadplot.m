function fH = sseegMakePrePostHeadplot()

%% Visualize

%% Choices to make:                                              


project_path   = '/Volumes/server/Projects/EEG/SSEEG/';
session_name   = 'SSEEG_20150403_wl_subj004';
session_prefix = 'Session_20150403_1145';

figureDir       = fullfile(project_path, 'Data', session_name, 'figures');
saveFigures     = true;     % Save figures in the figure folder?

% Load denoised data
load(fullfile(project_path, 'Data', session_name, 'processed', [session_prefix '_denoisedData_bb.mat']));
bbresults = results;

load(fullfile(project_path, 'Data', session_name, 'processed', [session_prefix '_denoisedData_sl.mat']));
slresults = results;

% Make a headplot of data before and denoising

figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(slresults.origmodel(1),icond, 'SNR');
    clims_sl = [0,25.6723];
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bbresults.origmodel(1),  icond, 'SNR');
    ab_snr2 = getsignalnoise(bbresults.finalmodel(1), icond, 'SNR');
    clims_ab = [0, max([ab_snr1, 12.4445])];
    

    % convert back into 157-channel space;
    newArray_absnr1 = nan(size(ab_snr1,1),128);
    newArray_absnr1(:,~badChannels) = ab_snr1;
    
    newArray_absnr2 = nan(size(ab_snr2,1),128);
    newArray_absnr2(:,~badChannels) = ab_snr2;
    
    newArray_slsnr1 = nan(size(sl_snr1,1),128);
    newArray_slsnr1(:,~badChannels) = sl_snr1;
    
       
    
    
    % plot spatial maps
    subplot(3,3,(icond-1)*3+1)
    plotOnEgi(newArray_slsnr1);
    set(gca,'CLim',[0 10]);
    colorbar;
    makeprettyaxes(gca,9,9);
    title(sprintf('SL no DN %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+2)
    plotOnEgi(newArray_absnr1);
    colorbar;

    set(gca,'CLim',[0 2]);
    makeprettyaxes(gca,9,9);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+3)
    plotOnEgi(newArray_absnr2);
    set(gca,'CLim',[0 2]);
    colorbar;
    makeprettyaxes(gca,9,9);
    title(sprintf('Broadband Post %s', condNames{icond}))
end

if saveFigures
    figurewrite(fullfile(figureDir,'figurePrePostDenoising'),[],0,'.',1);
end

