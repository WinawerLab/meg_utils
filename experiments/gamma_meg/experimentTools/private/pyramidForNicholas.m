%% Load an image
imageName = '/Volumes/server/Projects/MEG/Gamma/natural_images/nat_images_before/nat_image16.png';
imLoad = mean(double(imread(imageName)), 3);

mx = max(imLoad(:));
im = (imLoad - (mx/2))/mx;

%% Build a Laplacian pyramid
layers = 8;
[pyramid, LPR] = sketchBuildPyr(im, layers);

%% Thing 1: zero out the LPR entirely
LPR = zeros(size(LPR));

%% Another idea: zero-mean the LPR
% LPR = LPR - mean(LPR(:));
% LPR(LPR > 0.5) = 0.5; LPR(LPR < -0.5) = -0.5;

%% Thing 2: zero-mean most layers, pass out the lowest
% passCutoff = 7;
% for ii = 1:layers
%     if ii > passCutoff
%         pyramid{ii} = zeros(size(pyramid{ii}));
%     else
%         mn = mean(pyramid{ii}(:));
%         pyramid{ii} = pyramid{ii} - mn;
% 
%         % Clip overhang
%         pyramid{ii}(pyramid{ii} > 0.5) = 0.5;
%         pyramid{ii}(pyramid{ii} < -0.5) = -0.5;
%     end
% end

%% Rebuild the Laplacian pyramid
reconIm = sketchReconPyr(pyramid, LPR);
reconIm = reconIm - mean(reconIm(:));
reconIm(reconIm > 0.5) = 0.5;
reconIm(reconIm < -0.5) = -0.5;

%% Compare to zero-meaning the whole image at the pixel level
zeroMean = im - mean(im(:));
zeroMean(zeroMean > 0.5) = 0.5;
zeroMean(zeroMean < -0.5) = -0.5;

%% Show images
imRange = [-0.5, 0.5];
figure; clf; imshow(reconIm, imRange); title('Recon Image');
% figure(2); clf; imshow(zeroMean, imRange)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTRAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Desired new contrast settings
reconContrast = std(reconIm(:));
contrastFactor = 0.5;
targetContrast = reconContrast * contrastFactor;

%% Fiddle with the pixel-level image directly
basicConIm = reconIm * contrastFactor;

%% Fiddle with contrast across the pyramid

% Adapted recon script:
layers = length(pyramid);
recons = cell(1, layers);

lowPass = LPR;
for ii = layers:-1:1
    expandLP = imresize(lowPass, 2);
    lowPass = expandLP + pyramid{ii};
    
    % Fuss with the contrast;
    lowPass = lowPass - mean(lowPass(:));
    lowPass = (lowPass / (std(lowPass(:)))) * targetContrast;
    lowPass(lowPass > 0.5) = 0.5;
    lowPass(lowPass < -0.5) = -0.5;
    recons{ii} = lowPass;
end
pyrConIm = lowPass;

%% Comparison
figure(1); clf;
subplot(2, 3, 1); imshow(reconIm, imRange); title('Before contrast adjustment');
subplot(2, 3, 2); imshow(basicConIm, imRange); title('Basic contrast');
subplot(2, 3, 3); imshow(pyrConIm, imRange); title('Pyramid contrast');

subplot(2, 3, 4); hist(reconIm(:)); title(['STD: ', num2str(std(reconIm(:)))]); xlim(imRange);
subplot(2, 3, 5); hist(basicConIm(:)); title(['STD: ', num2str(std(basicConIm(:)))]); xlim(imRange);
subplot(2, 3, 6); hist(pyrConIm(:)); title(['STD: ', num2str(std(pyrConIm(:)))]); xlim(imRange);