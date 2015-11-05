function outputImage = normalizeIntensity(im)

%% import toolbox and images

% path to SHINE toolbox
projectDir = '~/matlab/SHINEtoolbox/'; 
addpath(genpath(projectDir))

sz = [768 768]; %desired pixel dimensions
maskRange = ones(sz); % no masking

visualizeBeforeAfter = true; 
useAverageStimulus = true; % use mean of 1000ms stims for template
% if false, uses uniform distribution


% image with uniform distribution of pixel intensity between 0:255
im_template = reshape(repmat(0:255, 2304,1), 768, 768);
im_template = sort(im_template(:)); % match function uses sorted list

if useAverageStimulus
    x = load('/Volumes/server/Projects/MEG/Gamma/stimuli/1000ms stimuli/gammaStimuli_params1.mat');
    stims = x.stimulus.images;
    im_template = mean(stims,3);
    im_template = sort(im_template);
end

% use saved images if not passed as arg
if nargin==0
    load('/Volumes/server/Projects/MEG/Gamma/stimuli/example_V1_electrode_faceshouses.mat');
    im = out.image;
    clear out;
    sz = 768;
end



scale   = sz/size(im,1); %rescale to MEG display res
toMatch = imresize(im,scale);
outputImage = ones(size(toMatch));

%% normalize

for i = 1:10
    
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
        
        pause(1);
        
    end
    
end

end

