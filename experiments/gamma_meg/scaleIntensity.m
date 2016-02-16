function outIm = scaleIntensity(im, targetContrast)
%% scaleIntensity
% scale the pixel values of a grayscale image to a specified standard 
% deviation given the pixel range -0.5 - 0.5
% @params: im - input image with px in the range [0 255]
%          targetContrast - desired SD, use > 1 for binarized
% @output: outIm - scaled image
% for use in MEG gamma experiment for gratings and noise patterns
% natural scene images use the function normalizeLaplacian 
% nicholas chua 2016

%% scale range
mx = 255; % assume 8 bits
im = (im - (mx/2))/mx; % scale range to [-0.5 0.5]

%% change contrast
equalizedImage = im/std(im(:)) * targetContrast;

% crop values that exceed the range
if targetContrast > 1,
    equalizedImage(equalizedImage>=0) =  0.5;
    equalizedImage(equalizedImage<0)  = -0.5;
else
    equalizedImage(equalizedImage>0.5)  = 0.5;
    equalizedImage(equalizedImage<-0.5) = -0.5;
end

%% scale back to [0 255] and cast to uint8
equalizedImage = (equalizedImage + 0.5)*mx;
outIm = uint8(equalizedImage);
end
