function [outIm]=normalizeLaplacian(im, varargin)
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
    imGray = double(imread('cameraman.tif')); % load image for demo
    targetContrast = 0.7;
else % use input image and contrast parameter
    imGray = mean(double(im),3);
    if isnumeric(varargin{1})
        targetContrast = varargin{1};
    else
        switch varargin{1}
            % change target contrast here
            case 'high'
                targetContrast = 5;%0.7;
            case 'medium'
                targetContrast = 0.1;
            case 'low'
                targetContrast = 0.02;
        end
    end
    
    if length(varargin)>1, mask = varargin{2}; end
end
    
if ~exist('mask', 'var'), mask = true(size(imGray)); end

mx = 255; % assume 8 bits
im = (imGray - (mx/2))/mx; % transform px intensities to range [-.5 .5]
imRange = [-.5 .5];

%% Build a Laplacian pyramid
layers = 8;
[origPyramid, LPR] = sketchBuildPyr(im, layers);
% pyramid contains the details lost at each level
% LPR is the low contrast residual 

if visualizeTransformations
   
    figure;
    for ii = 1:layers
        subplot(3,3,ii)
        imagesc(origPyramid{ii}, imRange); 
        colormap gray; axis image; title(sprintf('%ith pyramid',ii));        
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
reconIm = reconIm - median(reconIm(mask(:)));
equalizedImage = reconIm/std(reconIm(:)) * targetContrast;
%clip once more

if targetContrast > 1,
    equalizedImage(equalizedImage>=0) =  0.5;
    equalizedImage(equalizedImage<0)  = -0.5;
else
    equalizedImage(equalizedImage>0.5)  = 0.5;
    equalizedImage(equalizedImage<-0.5) = -0.5;
end

if visualizeTransformations
    figure(111)
    subplot(1,3,1);imshow(im, imRange); title('original image');
    subplot(1,3,2);imshow(reconIm, imRange); title('after high pass filter');
    subplot(1,3,3);imshow(equalizedImage, imRange); title('after equalization');
end

equalizedImage = (equalizedImage + 0.5)*mx; % rescale to original px range
outIm = uint8(equalizedImage);


end
