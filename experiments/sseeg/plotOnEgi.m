function varargout = plotOnEgi(data, showNumbers)
%plotOnEgi - Plots data on a standarized EGI net mesh
%function meshHandle = plotOnEgi(data)
%
%This function will plot data on the standardized EGI mesh with the
%arizona colormap.
%
%Data must be a 128 dimensional vector, but can have singleton dimensions,
%
%

if ~exist('showNumbers', 'var'), showNumbers = false; end

data = squeeze(data);
datSz = size(data);

if datSz(1)<datSz(2)
    data = data';
end

% Check for nans, if so: take median of all electrodes
for ii = 1:length(data)
    if isnan(data(ii))
        data(ii) = nanmedian(data);
    end
end

if size(data,1) == 128
tEpos = load('eeg128xy.mat');
tEpos = [ tEpos.xy, zeros(128,1) ];

tEGIfaces = mrC_EGInetFaces( false );

nChan = 128;
elseif size(data,1) == 256
    
tEpos = load('defaultFlatNet256.mat');
tEpos = [ tEpos.xy, zeros(256,1) ];

tEGIfaces = mrC_EGInetFaces256( false );
nChan = 256;
else
    error('Only good for 2 montages: Must input a 128 or 256 vector')
end


patchList = findobj(gca,'type','patch');
netList   = findobj(patchList,'UserData','plotOnEgi');

if isempty(netList),    
    handle = patch( 'Vertices', [ tEpos(1:nChan,1:2), zeros(nChan,1) ], ...
        'Faces', tEGIfaces,'EdgeColor', [0.5 0.5 0.5],  ...
        'FaceColor', 'interp');
    axis equal;
    axis off;
else
    handle = netList;
end

set(handle,'facevertexCdata',data,'linewidth',1,'markersize',20,'marker','.');
set(handle,'userdata','plotOnEgi');

% colormap(jmaColors('usadarkblue'));

if showNumbers,
    for ii = 1:size(tEpos,1); text(tEpos(ii,1), tEpos(ii,2), num2str(ii), 'FontSize', 10); end
end


if nargout >= 1
varargout{1} = handle;
end
