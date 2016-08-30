function [] = s_GAMMA_across_subjects_plot(sessions, denoised)
% creates figures out of the spectral data obtained from
% gamma_spectral_analysis.m
%
% gamma_visualize_meshes_across_subjects(sessions)
%
% Input:
%       sessions  :     Which sessions to plot?
%       deniosed  :     Flag to request denoised or undenoised data



%% Add fieldtrip path if necessary
if isempty(which('ft_prepare_layout'))
    addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'))
end

conditionNames = gamma_get_condition_names(sessions(1));
numConds       = length(conditionNames);
numChannels    = 157; % see if we can derive this from data instead of hardcode
summaryStat   = @(x) nanmean(x,3) ./ nanstd(x, [], 3);

% Get difference identity vectors (each condition - baseline)
contrasts = eye(numConds,numConds);

% Make baseline column negative to subtract from each cond
contrasts(:,find(strcmp(conditionNames,'Blank'))) = -1;

contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
numContrasts = size(contrasts,1);

% Predefine arrays
snrGamma = zeros(numContrasts, numChannels, length(sessions));
snrBroadband = zeros(numContrasts, numChannels, length(sessions));

% For plotting
scrsz = get(0,'ScreenSize');
threshold = 0;

for sessionNum = sessions
    
    sessionPath = meg_gamma_get_path(sessionNum);
    if denoised
        d = dir(fullfile(sessionPath,'processed','*_results*denoised*'));
        load(fullfile(sessionPath,'processed',d.name));
    else
        d = dir(fullfile(sessionPath,'processed','*_results*'));
        load(fullfile(sessionPath,'processed',d(1).name));
    end
    
    % Get opts
    opt = results.opt;
    opt.saveFigures = 1 ;
    opt.sessionPath = sessionPath; % Redefine in case data was analyzed on the HPC
    conditionNames = gamma_get_condition_names(opt.params.sessionNumber);
    
    % Define the file name of saved figures
    if opt.saveFigures
        if opt.MEGDenoise; postFix = '_denoised'; else postFix = ''; end
%         preFix = sprintf('s_%3d', opt.params.sessionNumber);
        saveName = sprintf('boots%d_%s_%s', opt.nBoot, postFix,...
            datestr(now, 'mm_dd_yy'));
    end
    
    %% Define summary metric (SNR or just mean)
    
    
    % Compute SNR
    tmpData = permute(results.broadbandPower, [2 1 3]);
    tmpData = reshape(tmpData, numConds, []);
    tmp = contrasts*tmpData;
    tmp = reshape(tmp, numContrasts, numChannels, opt.nBoot);
    snrBB = summaryStat(tmp)';
    
    tmpData = permute(results.gammaPower, [2 1 3]);
    tmpData = reshape(tmpData, numConds, []);
    tmp = contrasts*tmpData;
    tmp = reshape(tmp, numContrasts, numChannels, opt.nBoot);
    snrG  = summaryStat(tmp)';
    
    if length(snrBB) < 157; snrBB = to157chan(snrBB', ~opt.params.badChannels,0); end
    if length(snrG) < 157; snrG = to157chan(snrG', ~opt.params.badChannels,0); end
    
    % Add to the rest of the subjects
    snrBroadband(:,:,find(sessionNum==sessions)) = snrBB';
    
    snrGamma(:,:,find(sessionNum==sessions)) = snrG';
    
end



%% Plot meshes Gamma


for dataType = 1:2
    
    switch dataType
        case 1
            data = mean(snrGamma(1:12,:,:),3);
            str = 'Gamma';
        case 2
            data = mean(snrBroadband(1:12,:,:),3);
            str = 'Broadband';
    end
    
    figure; clf; set(gcf, 'position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)]);  set(gcf, 'name', sprintf('%s SNR', str))
    plot_range = [-1 1] * ceil(max(abs(data(:))));
    
    for c = 1:12
        subplot(4,3,c)
        dataToPlot = data(c,:);
        dataToPlot(abs(dataToPlot) < threshold) = 0;
        ft_plotOnMesh(dataToPlot, conditionNames{c});
        set(gca, 'CLim', plot_range);
        colormap(jmaColors('coolhotcortex'))
    end
    
    if opt.saveFigures;
        if ~exist(fullfile(fileparts(opt.sessionPath),'AcrossSubjects','meshes'),'dir'); mkdir(fullfile(fileparts(opt.sessionPath),'AcrossSubjects','meshes')); end
        hgexport(gcf, fullfile(fileparts(opt.sessionPath),'AcrossSubjects','meshes',sprintf('SNR_mesh_%s_%s_%d.eps', str, saveName, threshold)));
    end
    
end





