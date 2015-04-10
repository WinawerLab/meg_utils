function fH = sseeg_plotOnMesh(sensor_data, title_txt, figure_num, plotType)
%Plot a surface map of the MEG sensor data
% fH = ssm_plotOnMesh(sensor_data, [title_txt], [figure_num], meg_files, [plotType])
%
% Inputs
%   sensor_data: 1x128 vector of sensor data (if longer than 157, then
%           truncate)
%   title: optional string for title of plot
%   figure_num: optional integer to specify figure number
%   plotType: any of '2d', '3d', 'both'
%
% Output
%   fH: figure handle (integer)
%
% Example:
%   fH = sseeg_plotOnMesh(randn(1,128), 'my random data', 999);


% check inputs
if ~exist('plotType', 'var') || isempty(plotType), plotType = 'both'; end

% check length of sensor data
if length(sensor_data) > 128, sensor_data = sensor_data(1:28); end



% Load the sensor positions in  space 
tmp = load('meg160xyz');
xyz = tmp.xyz;


%% 3D plot - Requires fieldtrip

% add fieldtrip matlab code
if isempty(which('ft_analysispipeline')),
    addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'));
end

switch lower(plotType)
    case {'3d', 'both'}
        % set up figure
        if exist('figure_num', 'var'), fH = figure(figure_num); clf;
        else                           fH = figure; clf; end
        set(fH, 'Color', 'w')
        
        % Get labels
        colorbar
        ylabel('    Right       Left     ')
        xlabel('    Posterior       Anterior     ')
        zlabel('    Inferior       Superior     ')
        
        ft_plot_topo3d(xyz,sensor_data); hold on;
        label_add(xyz)
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
end

%% 2D TOPOPLOT - Requires fieldtrip

switch lower(plotType)
    case {'2d', 'both'}
        
        
%         data_hdr = ft_read_header(fullfile('raw',meg_files(1).name));
        
        if notDefined('data_hdr')
           data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
        end
        
        cfg=[];
        
        cfg.layout          = ft_prepare_layout(cfg, data_hdr);
        cfg.style           ='straight';
        % cfg.style         ='blank';
        cfg.electrodes      ='off';
        cfg.colorbar        ='no';
        cfg.maplimits       ='maxmin';
        cfg.data            = sensor_data';
        cfg.interpolation   = 'v4';


        
        for ii = 1:length(cfg.data)
            if isnan(cfg.data(ii))
                cfg.data(ii) = nanmedian(sensor_data);
            end
        end
        figure; clf;
        topoplot(cfg,cfg.data)
        
        fH = gcf;
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
        
end

