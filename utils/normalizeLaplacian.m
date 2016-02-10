function [outIm]=normalizeLaplacian(im, args)
%% normalizeLaplacian
% normalize the contrast level of an array of images
% by 1)high pass filtering of using Laplacian Pyramid
%    2)equalize the standard deviation of px intensity to a certain amount
% inputs: im - 3D array (sz x sz x n) of grey scale images
%         arg - string ('high' 'medium' or 'low')
%%%%%%%%%%% Leave arguments blank for a demonstration %%%%%%%%%%%%%%%%%%%%%

%% parameters and arguments

visualizeTransformations = false;

if nargin == 0 % load example image
    imageName = '/Volumes/server/Projects/MEG/Gamma/natural_images_tools/nat_images_before/nat_image16.png';
    imLoad = mean(double(imread(imageName)), 3); % load image for demo
    targetContrast = 0.7;
else % use input image and contrast parameter
    imLoad = mean(double(im),3);
    switch args
        % change target contrast here
        case 'high'
            targetContrast = 0.7;
        case 'medium'
            targetContrast = 0.1;
        case 'low'
            targetContrast = 0.02;
    end
end
    
mx = max(imLoad(:)); 
im = (imLoad - (mx/2))/mx; % transform px intensities to range [-.5 .5]
imRange = [-.5 .5];

%% add path to functions
% sketchReconPyr.m and sketchBuildPyr.m
path = '//Volumes/server/Projects/MEG/Gamma/natural_images_tools/pyramid';
addpath(genpath(path));

%% Build a Laplacian pyramid
layers = 8;
[origPyramid, LPR] = sketchBuildPyr(im, layers);
% pyramid contains the details lost at each level
% LPR is the low contrast residual 

if visualizeTransformations
   
    for ii = 1:layers
        figure(3);imagesc(origPyramid{ii}, imRange); colormap gray; axis image; title(sprintf('%ith pyramid',ii));
        pause(0.5);
    end
      
end

%% Reconstruct with zero-ed out low pass residual

% zero out LPR
LPR = zeros(size(LPR));

% Reconstruct image from pyramidLayers and LPR
reconIm = sketchReconPyr(origPyramid, LPR);
reconIm = reconIm - mean(reconIm(:)); % zero-mean again
reconIm(reconIm > 0.5) = 0.5;
reconIm(reconIm < -0.5) = -0.5;


%% Change contrast

equalizedImage = reconIm/std(reconIm(:)) * targetContrast;
%clip once more
equalizedImage(equalizedImage>0.5) = 0.5;
equalizedImage(equalizedImage<-0.5) = -0.5;

if visualizeTransformations
    figure(111)
    subplot(1,3,1);imshow(im, imRange); title('original image');
    subplot(1,3,2);imshow(reconIm, imRange); title('after high pass filter');
    subplot(1,3,3);imshow(equalizedImage, imRange); title('after equalization');
end

equalizedImage = (equalizedImage + 0.5)*mx; % rescale to original px range
outIm = uint8(equalizedImage);


end
