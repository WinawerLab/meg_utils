
%% Normalize contrast of images using Laplacian Pyramid
% work in progress
% Cathetine Olsson & Nicholas Chua

% TODO:
% - implement as function


%% Load images (remove and turn into input parameter)
imageName = '/Volumes/server/Projects/MEG/Gamma/natural_images/nat_images_before/nat_image16.png';
imLoad = mean(double(imread(imageName)), 3);

mx = max(imLoad(:));
im = (imLoad - (mx/2))/mx; %im has pixel range between -.5 and 5

%% Parameters
visualizeTransformations = true;
targetStandardDev = 0.1;


%%
%%%%%%%%%%%%%%%%%%% Remove Low Spatial Freq Detail %%%%%%%%%%%%%%%%%%%%%
%% Use high pass filter

imUse = im - mean(im(:)); % make the mean pixel intensity 0
if mod(size(imUse, 1), 2) == 0
    imUse = imresize(im, (size(imUse, 1)-1)/size(imUse,1));
end

sz = size(imUse);

% FFT of image
fftIm = fftshift(fft2(imUse));
logAbsFft = log(abs(fftIm));

% Build a filter / mask
mid = ceil((sz(1)+1)/2);
lin = linspace(-1,1,sz(1));
[x, y] = meshgrid(lin,lin);
r = sqrt(x.^2 +y.^2);

mask = ones(sz);
mask(r<0.04) = 0;

% blur mask to remove edge artifacts from 2D fft
gausBlur = fspecial('gaussian', sz,15);
blurredMask = imfilter(mask, gausBlur, 'replicate');

% Apply the filter
maskedFFT = abs(blurredMask.*fftIm);

% Reconstruct the image
highPassedIm = ifft2(ifftshift(maskedFFT));



%% Visualize 

imRange = [-0.5, 0.5];
figure(1); 
subplot(2,1,1); imshow(imUse, imRange); title('before HPF');
subplot(2,1,2); imshow(highPassedIm, imRange); title('High Pass Filtered');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Contrast %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Equalize variance to target

% change std to target stf
equalizedIm = highPassedIm/std(highPassedIm(:)) * targetStandardDev;

% clip to range
equalizedIm(equalizedIm > 0.5) = 0.5;
equalizedIm(equalizedIm < -0.5) = -0.5;

imRange = [-0.5, 0.5];
figure(2); 
subplot(2,2,1); imshow(highPassedIm, imRange); title('before equalization');
subplot(2,2,2); imshow(equalizedIm, imRange); title('after equalization');
subplot(2,2,3); hist(highPassedIm(:)); title(sprintf('STD: %d', std(highPassedIm))); xlim(imRange);
subplot(2,2,4); hist(abs(equalizedIm(:))); title(sprintf('STD: %d', std(equalizedIm)));xlim(imRange);
