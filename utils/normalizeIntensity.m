function outputImage = normalizeIntensity(im)
% normalizes the intensity of an array of images based on the px intensity
% of a template image

%% import toolbox and images

% path to SHINE toolbox
projectDir = '~/matlab/SHINEtoolbox/'; 
addpath(genpath(projectDir))

% which image to use as template for intensity?
useAverageStimulus = true; % use mean of 1000ms stims for template
% if false, uses uniform distribution

visualizeBeforeAfter = true; 

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
end

% match function requires a mask
maskRange = ones(sz); % no masking

if useAverageStimulus % use average intensities of the previous experiment's stimuli
    x = load('/Volumes/server/Projects/MEG/Gamma/stimuli/1000ms stimuli/gammaStimuli_params1.mat');
    stims = x.stimulus.images;
    im_template = mean(stims,3);
    im_template = sort(im_template);
else
    % image with uniform distribution of pixel intensity between 0:255
    im_template = reshape(repmat(0:255, 2304,1), 768, 768);
    im_template = sort(im_template(:)); % match function uses sorted list
end


 %rescale to MEG display res

outputImage = ones(size(toMatch));

%% normalize

for i = 1:size(im,3)
    
    thisImage = toMatch(:,:,i);
    
    outputImage(:,:,i) = match(thisImage, im_template, maskRange);
    
end

%% visualize

if visualizeBeforeAfter
    
    for ii = 1:size(outputImage, 3)
        
        figure (2), clf, colormap gray
        
        subplot(2, 2, 1)
        imagesc(maskRange)
        
        subplot(2, 2, 2)
        imagesc(toMatch(:,:,ii))
        
        subplot(2, 2, 3)
        imagesc(im_template)
        
        subplot(2, 2, 4)
        imagesc(outputImage(:,:,ii))
        
        waitforbuttonpress;
        
    end
    
end

end

