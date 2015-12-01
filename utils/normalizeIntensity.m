
% normalizes the intensity of an array of images (sz x sz x N) based on the px intensity
% of a template image

% To adjust (could be optional inputs)

% which image to use as template for intensity? if false, uses uniform distribution
useAverageStimulus = false; % true = use mean of 1000ms stims for template

visualizeBeforeAfter = true; 

normalizeMethod = 'rms'; % rms or histogram

low_or_high_contrast = 'low'; 

masking = true;

save_images = false; 
%% import toolbox and images

% path to SHINE toolbox
projectDir = '~/matlab/SHINEtoolbox/'; 
addpath(genpath(projectDir))



% use saved images if not passed as arg
if nargin==0 % load all images from Hermes et al. 
    load('/Volumes/server/Projects/MEG/Gamma/stimuli/example_V1_electrode_faceshouses.mat');
    im = out.image;
    clear out;
    sz = 768; % scale to MEG display
    scale   = sz/size(im,1);
    toMatch = imresize(im,scale);
else % use input images
    sz = size(im,3); 
    toMatch = im;
end


switch normalizeMethod
    case 'histogram'
        if useAverageStimulus % use average intensities of the previous experiment's stimuli
            x = load('/Volumes/server/Projects/MEG/Gamma/stimuli/1000ms stimuli/gammaStimuli_params1.mat');
            stims = x.stimulus.images;
            im_template = mean(stims,3);
            im_template = sort(im_template);
        else
            % image with uniform distribution of pixel intensity between 0:255
            switch low_or_high_contrast
                case 'high'
                    im_template = reshape(repmat(0:255, 2304,1), 768, 768);
                case 'medium'
                    im_template = reshape(repmat(round(linspace(128-64, 128+64, 256)), 2304,1), 768, 768);
                case 'low'
                    im_template = reshape(repmat(round(linspace(128-28, 128+28, 256)), 2304,1), 768, 768);
            end
            im_template = sort(im_template(:)); % match function uses sorted list
        end
        
    case 'rms'
        switch low_or_high_contrast
            case 'high',    im_template = 128 + 1000*[-1 1];
            case 'medium',  im_template = 128 + 64*[-1 1];
            case 'low',     im_template = 128 + 14*[-1 1];
        end
        template_sd = std(im_template);
end
 %rescale to MEG display res

outputImage = ones(size(toMatch));

%% masking

% if masking
%     v           = linspace(-1, 1, 768);
%     [x, y]      = meshgrid(v, v);
%     mask = sqrt(0.2*x.^2 + 0.2*y.^2)< 0.6;
%     for i = 1:size(toMatch,3)
%         toMatch(:,:,i) = toMatch(:,:,i) .* mask;
%     end
% end






%% normalize
fprintf('[%s]: Normalizing %d images ', mfilename, size(im,3)); 
for i = 1:size(im,3)
    thisImage = toMatch(:,:,i);
    
    switch normalizeMethod
        case 'histogram'            
            outputImage(:,:,i) = match(thisImage, im_template, mask);
        case 'rms'            
            outputImage(:,:,i) = thisImage/std(thisImage(:))*template_sd;
    end
    fprintf('.'); drawnow
end

fprintf('\n')


%% Make a mask 
% 
[x, y] = meshgrid(linspace(-1,1,sz));
% mask = ones(sz); % no masking

R             = sqrt(0.2*x.^2 + 0.2*y.^2);
r_min         = .1;
r_max         = .5;
Edge          = (R-r_min) / (r_max - r_min);
Edge(R<r_min) = 0;
Edge(R>r_max) = 1;

cosMask       = (cos(Edge*pi)+1)/2;
blank         = ones(sz,sz) * 0.5;


maskedOutputImage(:,:,1) = outputImage(:,:,1).*cosMask + (1-cosMask).*blank;

% understand cosine mask
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


subplot(3, 2, 4)
imagesc(blank)
title('blank')

subplot(3, 2,5)
imagesc(outputImage(:,:,1))
title('beforemask')

subplot(3, 2,6)
imagesc(maskedOutputImage(:,:,1))
title('aftermask')

%%

% saveLocText  = 'Images/ellipseMaskedFaceStimulus';
%mkdir(saveLocText)
% saveLocation = sprintf('%s%s', project_pth, saveLocText);
% cd(saveLocation) 
% delete('*.png')

% 
% for kk = 1 : size(outputImage,3)
%     maskedOutputImage(:,:,kk) = outputImage(:,:,kk).*cosMask + (1-cosMask).*blank;
%     figure (3), colormap gray
%     subplot(1, 3, 1), imagesc(outputImage(:,:,kk)); axis square
%     subplot(1, 3, 2), imagesc(maskedOutputImage(:,:,kk)); axis square
%     subplot(1, 3, 3), imagesc(maskedOutputImage(:,:,kk)); axis square
%   %  imwrite(uint8(imresize(double(maskedFace{kk}),[256, 256])), ...
%   %      sprintf('ellipseMaskedface%d.png',  kk));
%     pause (1)
% end
% 


%% save/visualize

outputImage = outputImage + 128;
outputImage = uint8(outputImage);
toMatch     = toMatch + 128;
toMatch     = uint8(toMatch);

if save_images
for ii = 1:size(outputImage,3)
    file_name = sprintf('%s_contrast_image number_%d.png', low_or_high_contrast, ii);
    test_save_pth = '/Volumes/server/Projects/MEG/Gamma/natural_images/nat_images_rms/circular_masked/';
    imwrite(outputImage(:,:,ii), strcat(test_save_pth, file_name));
end
end


if visualizeBeforeAfter
    
    % Visualize after
    figure,
    a = reshape(reshape(1:72,12,6)', 1, 72);
    for jj = 1:size(outputImage, 3),
        subplot(6,12,a(jj)), imagesc(uint8(outputImage(:,:,jj)), [0 255]),
        axis image off
        colormap gray
    end
    

    
    for ii = 1:size(outputImage, 3)
        
            

        figure (2), clf, colormap gray
        set(gcf, 'Name', sprintf('Image %d', ii));
        subplot(2, 2, 1)        
        imagesc(cosMask, [0 1]), title('Mask')
        
        subplot(2, 2, 2)
        imagesc(toMatch(:,:,ii)), title('Original')
        
        subplot(2, 2, 3)
        imagesc(im_template, [0 255]), title('Template')
        
        subplot(2, 2, 4)
        imagesc(outputImage(:,:,ii), [0 255]), title('Normalized')
        
        waitforbuttonpress;
        
    end
    
end



