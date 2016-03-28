function [maskedIm, isolatedStd, isolatedMn] = cosineMask(im, r_min,ii)
%% cosineMask
% applies a soft edge cosine mask to a square 2D image
% change visualize to true to learn about the mask

if ~exist('r_min', 'var'), r_min = .8; end

visualize = false;
saveimage = true;

projectPath = '/Volumes/server/Projects/MEG/Gamma';
savePath    = fullfile(projectPath, 'stimuli/natural_images/histograms3.22');
%% make the mask
sz = size(im,1);
bg = 128; % the mask should be grey (mean luminanace)

[x, y]        = meshgrid(linspace(-1,1,sz));
R             = sqrt(x.^2 + y.^2);
Edge          = (R-r_min) / (1 - r_min);
Edge(R<r_min) = 0;
Edge(R>1)     = 1;

tmp = double(im)-bg;

cosMask = (cos(Edge*pi)+1)/2 ;

%% apply mask
maskedIm = uint8(tmp .*cosMask + bg); 

%% measure RMS inside mask
isolatedImage = double(maskedIm);
isolatedImage(R > r_min) = NaN;

isolatedStd = nanstd(isolatedImage(:));
isolatedMn  = nanmean(isolatedImage(:));

%% save just image
if saveimage && mod(ii,9) < 4
    figure (100), clf, colormap gray
    imshow(maskedIm)
    hgexport(gcf,fullfile(savePath, sprintf('/image_%d', ii)));
end

%% Visualization
if visualize
%if ii > 27 && ii < 58
    figure (100), clf, colormap gray
%     subplot(3, 2, 1)
%     imagesc(R),
%     title('R')
%     
%     subplot(3, 2, 2)
%     imagesc(Edge)
%     title('Edge')
%     
%     subplot(3, 2, 3)
%     imagesc(cosMask)
%     title('cosMask')
%     
%     subplot(3, 2,4)
%     imagesc(im, [1 225])
%     title('beforemask')
    
    subplot(1,2,1)% subplot(3, 2,5)
    imshow(maskedIm)
    title('aftermask')
    
    subplot(1,2,2)% subplot(3,2,6)
    hist(isolatedImage(:))
    title(sprintf('distribution SD = %2.2f MN = %2.2f', isolatedStd, isolatedMn))
    
    hgexport(gcf,fullfile(savePath, sprintf('/histogram_%d', ii)));
    %waitforbuttonpress;
    pause(0.1);
%end
end
end