function fH = ft_plotOnMesh(sensor_data, title_txt, figure_num, plotType, varargin)
%Plot a surface map of the MEG sensor data
% fH = ft_plotOnMesh(sensor_data, [title_txt], [figure_num], [plotType])
%
% Inputs
%   sensor_data: 1x157 vector of sensor data (if longer than 157, then
%           truncate)
%   title_txt: optional string for title of plot
%   figure_num: optional integer to specify figure number
%   plotType: '2d' or '3d'
%
% Output
%   fH: figure handle (integer)
%
% Example:
%   fH = ft_plotOnMesh(randn(1,157), 'my random data', 999,'2d');


% check inputs
if notDefined('plotType'), plotType = '2d'; end

% check length of sensor data
if length(sensor_data) > 157, sensor_data = sensor_data(1:157); end

% set up figure
if notDefined('figure_num'), fH = gcf;
else fH = figure(figure_num); clf; set(fH, 'Color', 'w'), end



% % add fieldtrip matlab code
% if isempty(which('ft_analysispipeline')),
%     addpath(genpath('/Volumes/server/Projects/MEG/code/fieldtrip'));
% end

switch lower(plotType)
    case '3d'
        
        % 3D plot - Requires fieldtrip
        
        % Load the sensor positions in 3 space
        tmp = load('meg160xyz');
        xyz = tmp.xyz;
        
        % Get labels
        colorbar
        ylabel('    Right       Left     ')
        xlabel('    Posterior       Anterior     ')
        zlabel('    Inferior       Superior     ')
        
        ft_plot_topo3d(xyz,sensor_data); hold on;
        % label_add(xyz)
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
        
    case '2d'
        % 2D TOPOPLOT - Requires fieldtrip
        
        if ~exist('data_hdr', 'var')
            data_hdr = load('hdr'); data_hdr = data_hdr.hdr;
        end
        
        % Default CFG
        cfg=[];        
        cfg.layout          = ft_prepare_layout(cfg, data_hdr);
        cfg.style           ='straight';
        %cfg.electrodes      ='numbers';
        %cfg.electrodes      = 'dotnum';
        cfg.electrodes      = 'on';
        cfg.emarkersize     = 5;
        cfg.efontsize       = 15;
        cfg.ecolor          = 'k';
        cfg.colorbar        ='yes';
        cfg.maplimits       ='maxmin';
        cfg.interpolation   = 'v4';

        cfg.data            = sensor_data';
        
        % Add custom options if requested
        if exist('varargin', 'var') 
           for ii = 1:2:length(varargin)
              switch lower(varargin{ii}) 
                  case 'electrodes', cfg.electrodes = varargin{ii+1};
                  case 'interpolation', cfg.interpolation = varargin{ii+1};                      
                  case 'colorbar', cfg.colorbar = varargin{ii+1};
                  case 'clim', clim = varargin{ii+1};
              end    
           end            
        end
        
        for ii = 1:length(cfg.data)
            if isnan(cfg.data(ii))
                cfg.data(ii) = nanmedian(sensor_data);
            end
        end
        topoplot(cfg,cfg.data)
        
        % add a title if requested
        if exist('title_txt', 'var') && ~isempty(title_txt), title(title_txt); end
        if exist('clim', 'var'), set(gca, 'CLim', clim); end
end

