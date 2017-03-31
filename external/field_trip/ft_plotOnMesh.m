function fH = ft_plotOnMesh(sensorData, titleTxt, figureNum, plotType, varargin)
%Plot a surface map of the MEG sensor data
% fH = ft_plotOnMesh(sensor_data, [title_txt], [figure_num], [plotType])
%
% Inputs
%   sensorData: 1x157 vector of sensor data (% if longer or shorter, it
%   will check for other layouts, ie. Yogokawa Con or Neuromag 360
%   title_txt: optional string for title of plot
%   figure_num: optional integer to specify figure number
%   plotType: '2d' or '3d'
%
%   Requires the Fieldtrip toolbox
%
% Output
%   fH: figure handle (integer)
%
% Example:
%   fH = ft_plotOnMesh(randn(1,157), 'my random data', 999,'2d');


%% Check inputs

% Plot Type
if notDefined('plotType'), plotType = '2d'; end

% Length of sensor data and layout
cfg = [];

if length(sensorData) <= 102 % Combined NeuroMag 204 channel MEG
        data_hdr = load('neuromag360_sample_hdr_combined.mat'); data_hdr = data_hdr.hdr;
        cfg.layout = ft_prepare_layout(cfg, data_hdr); 
    elseif length(sensorData) <= 104 % Combined Yokogawa 208 channel MEG
        data_hdr = load('yokogawa_con_example_hdr.mat'); data_hdr = data_hdr.hdr;
        cfg.layout = ft_prepare_layout(cfg, data_hdr); % Create our own layout with Fieltrip function
    elseif length(sensorData) <= 157 % Uncombined/standard Yokogawa 157 channel MEG
        data_hdr = load('meg160_example_hdr.mat'); data_hdr = data_hdr.hdr;
        cfg.layout = ft_prepare_layout(cfg, data_hdr); % Create our own layout with Fieltrip function
    elseif length(sensorData) <= 204 % Uncombined NeuroMag 204 planar channel MEG
        data_hdr = load('neuromag360_sample_hdr.mat'); data_hdr = data_hdr.hdr;
        cfg.layout = ft_prepare_layout(cfg, data_hdr); 
    elseif length(sensorData) <= 208
        data_hdr = load('yokogawa_con_example_hdr.mat'); data_hdr = data_hdr.hdr;
        cfg.layout = ft_prepare_layout(cfg, data_hdr);
    else
        error(sprintf('(%s): Can''t find MEG layout configuration', mfilename))
end

% Set up figure
if notDefined('figureNum'), fH = gcf;
else fH = figure(figureNum); clf; set(fH, 'Color', 'w'), end

% Plot for each type (
switch lower(plotType)
    case '3d'
        
        % 3D plot - Requires fieldtrip
        
        % Load the sensor positions in 3 space
        xyz = data_hdr.grad.chanpos;
        
        % Get labels
        colorbar
        ylabel('    Right       Left     ')
        xlabel('    Posterior       Anterior     ')
        zlabel('    Inferior       Superior     ')
        
        ft_plot_topo3d(xyz,sensorData); hold on;
        
        % add a title if requested
        if exist('titleTxt', 'var') && ~isempty(titleTxt), title(titleTxt); end
        
    case '2d'
        % 2D TOPOPLOT - Requires fieldtrip     
        chanX  = cfg.layout.pos(1:length(sensorData),1);
        chanY  = cfg.layout.pos(1:length(sensorData),2);
        cfg.data = sensorData;
        
        opt = {'interpmethod','v4',... How to interpolate the data?
        'interplim','mask',... Mask the data such that it doesn't exceed the outline
        'gridscale',170,... How fine do you want the grid?
        'outline',cfg.layout.outline,... Create the lines of the head, nose and ears
        'shading','flat', ... How to interpolate the in the outline
        'mask',cfg.layout.mask,... 
        'datmask', []};
        
        if exist('varargin', 'var')
            for ii = 1:2:length(varargin)
                % paired parameter and value
                parname = varargin{ii};
                val     = varargin{ii+1};

                % check wehther this parameter exists in the defaults
                existingparnames = opt(1:2:end);
                tmp = strfind(existingparnames, parname);
                idx = find(cellfun(@(x) ~isempty(x), tmp));

                % if so, replace it; if not add it to the end of opt
                if ~isempty(idx), opt{idx+1} = val;
                else, opt{end+1} = parname; opt{end+1} = val; end
            end
        end
        
        % Deal with NaNs
        for ii = 1:length(cfg.data)
            if isnan(cfg.data(ii))
                cfg.data(ii) = nanmedian(sensor_data);
            end
        end
        
        % Plot it!
        ft_plot_topo(chanX,chanY,cfg.data, opt{:});

        % Make plot pretty
        ft_plot_lay(cfg.layout,...
            'box','no',...
            'label','no',...
            'point','yes', ...
            'pointsymbol','.',...
            'pointcolor','k',...
            'pointsize',8, ...
            'colorbar', 'yes', ...
            'style','straight', ...
            'maplimits', 'maxmin')

        % add a title if requested
        if exist('titleTxt', 'var') && ~isempty(titleTxt), title(titleTxt); end
        if exist('clim', 'var'), set(gca, 'CLim', clim); end
end

