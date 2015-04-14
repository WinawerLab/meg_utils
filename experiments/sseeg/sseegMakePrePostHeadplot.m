function fH = sseegMakePrePostHeadplot(projectPath,sessionName,sessionPrefix,saveFigures)

%% Visualize  Stimulus signal, broadband before and after denoising
% We take the bootstrapped beta values for each condition (in this case we
% defined 3 (Full, left, right) for one subject. We calculate the SNR by
% taking the median of the bootstrapped beta values and divide this by the
% standard error. We pad the SNR data with NaNs for the missing bad
% channels and then plot it on a flat topoplot of the head.

% INPUTS:
% projectPath      : String variable, path where saved denoised results are stored
%                       Assumes that there is a folder called
%                       'processed' with mat files containing the outputs
%                       of the denoising algorithm.
% sessionName      : Folder name of subject
% sessionPrefix    : Name of particular session.
% saveFigures       : Boolean to save figure or not

% OUTPUTS:
% fH                : Figure handle of headlplot

% Dependencies
% This function depends on meg_utils and denoiseproject repositories.

%% Options to set:                                              
if notDefined('projectPath'),   projectPath = '/Volumes/server/Projects/EEG/SSEEG/'; end
if notDefined('sessionName'),   sessionName = 'SSEEG_20150403_wl_subj004'; end
if notDefined('sessionPrefix'), sessionPrefix = 'Session_20150403_1145'; end
if notDefined('saveFigures'),   saveFigures = false; end

figureDir       = fullfile(projectPath, 'Data', sessionName, 'figures');
if ~exist(figureDir,'dir'); mkdir(figureDir); end

%% Load denoised data
load(fullfile(projectPath, 'Data', sessionName, 'processed', [sessionPrefix '_denoisedData_bb.mat']));
bbResults = results;

load(fullfile(projectPath, 'Data', sessionName, 'processed', [sessionPrefix '_denoisedData_sl.mat']));
slResults = results;

%% Make a headplot of data before and denoising
figure('position',[1,600,1400,800]);
condNames = {'Stim Full','Stim Left','Stim Right'};
for icond = 1:3
    % get stimulus-locked snr
    tmpSLSNR1 = getsignalnoise(slResults.origmodel(1),icond, 'SNR');
    climsSL = [0,10];
    % get broadband snr for before and after denoising
    tmpBBSNR1 = getsignalnoise(bbResults.origmodel(1),  icond, 'SNR');
    tmpBBSNR2 = getsignalnoise(bbResults.finalmodel(1), icond, 'SNR');
    climsBB = [0, max([tmpBBSNR1, tmpBBSNR2])];

    % convert back into 128-electrode space
    slsnr1 = nan(size(tmp_sl_snr1,1),128);
    slsnr1(:,~badChannels) = tmpSLSNR1;
    
    absnr1 = nan(size(tmp_ab_snr1,1),128);
    absnr1(:,~badChannels) = tmpBBSNR1;
    
    absnr2 = nan(size(tmp_ab_snr2,1),128);
    absnr2(:,~badChannels) = tmpBBSNR2;
    
    % plot spatial maps
    subplot(3,3,(icond-1)*3+1)
    plotOnEgi(slsnr1);
    set(gca,'CLim',climsSL);
    colorbar;
    makeprettyaxes(gca,9,9);
    title(sprintf('SL no DN %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+2)
    plotOnEgi(absnr1);
    colorbar;
    set(gca,'CLim',climsBB);
    makeprettyaxes(gca,9,9);
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(3,3,(icond-1)*3+3)
    plotOnEgi(absnr2);
    set(gca,'CLim',climsBB);
    colorbar;
    makeprettyaxes(gca,9,9);
    title(sprintf('Broadband Post %s', condNames{icond}))
end

if saveFigures
    figurewrite(fullfile(figureDir,'figurePrePostDenoising'),[],0,'.',1);
end

