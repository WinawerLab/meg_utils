function varargout = plotOnEgi(data)
%plotOnEgi - Plots data on a standarized EGI net mesh
%function meshHandle = plotOnEgi(data)
%
%This function will plot data on the standardized EGI mesh with the
%arizona colormap.
%
%Data must be a 128 dimensional vector, but can have singleton dimensions,
%
%
data = squeeze(data);
datSz = size(data);

if datSz(1)<datSz(2)
    data = data';
end

if size(data,1) == 128
tEpos = load('defaultFlatNet.mat');
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
        'Faces', tEGIfaces,'EdgeColor', [ 0.5 0.5 0.5 ], ...
        'FaceColor', 'interp');
    axis equal;
    axis off;
else
    handle = netList;
end

set(handle,'facevertexCdata',data,'linewidth',1,'markersize',20,'marker','.');
set(handle,'userdata','plotOnEgi');

colormap(jmaColors('usadarkblue'));

if nargout >= 1
varargout{1} = handle;
end
