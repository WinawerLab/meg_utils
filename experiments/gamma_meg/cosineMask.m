function maskedIm = cosineMask(im)
%% cosineMask
% applies a soft edge cosine mask to a square 2D image
% change visualize to true to learn about the mask

visualize = true;
%% make the mask
sz = size(im,1); 

[x, y]        = meshgrid(linspace(-1,1,sz));
R             = sqrt(x.^2 + y.^2);
r_min         = .6; % adjust this to change size of aperture
Edge          = (R-r_min) / (1 - r_min);
Edge(R<r_min) = 0;
Edge(R>1)     = 1;

cosMask = (cos(Edge*pi)+1)/2;

%% apply mask
maskedIm = double(im).*cosMask; % int will cause a round up
maskedIm = uint8(maskedIm); % cast back to uint8
%% Visualization
if visualize
    figure (100), clf, colormap gray
    subplot(3, 2, 1)
    imagesc(R),
    title('R')
    
    subplot(3, 2, 2)
    imagesc(Edge)
    title('Edge')
    
    subplot(3, 2, 3)
    imagesc(cosMask)
    title('cosMask')
    
    subplot(3, 2,5)
    imagesc(im, [1 225])
    title('beforemask')
    
    subplot(3, 2,6)
    imagesc(maskedIm, [1 225])
    title('aftermask')
end
end