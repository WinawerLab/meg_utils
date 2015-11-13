function outputImage = normalizeIntensity(im)
% normalizes the intensity of an array of images based on the px intensity
% of a template image

% To adjust (could be optional inputs)

% which image to use as template for intensity? if false, uses uniform distribution
useAverageStimulus = false; % use mean of 1000ms stims for template

visualizeBeforeAfter = true; 

normalizeMethod = 'rms'; % rms or histogram

low_or_high_contrast = 'low'; 
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

% match function requires a mask
[x, y] = meshgrid(linspace(-1,1,sz));
% mask = ones(sz); % no masking
mask = x.^2 + y.^2 < 1;

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
%% visualize

outputImage = outputImage + 128;
outputImage = uint8(outputImage);
toMatch     = toMatch + 128;
toMatch     = uint8(toMatch);

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
        imagesc(mask, [0 1]), title('Mask')
        
        subplot(2, 2, 2)
        imagesc(toMatch(:,:,ii)), title('Original')
        
        subplot(2, 2, 3)
        imagesc(im_template, [0 255]), title('Template')
        
        subplot(2, 2, 4)
        imagesc(outputImage(:,:,ii), [0 255]), title('Normalized')
        
        waitforbuttonpress;
        
    end
    
end

end

