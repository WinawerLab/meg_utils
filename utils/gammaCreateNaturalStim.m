% Code to create gamma stimuli
% includes natural images, gratings,
% and high contrast noise

%% %% Stimulus parameters

bckgrnd      = 128;      % this will be the background and mean pixel value
                        % therefore our images should range from 1 to 255
rng          = [1 255];  % range of image values
image_dur    = 1.0;      % seconds
Blank_dur    = 0.5;

scale_images = @(x) uint8((x - min(x(:))) / (max(x(:)) - min(x(:))) * diff(rng) + min(rng));
nImages      = 7; % images per category
nCategories  = 10;
nNatural     = 7; % number of natural images

nTotal       = nImages * nCategories;
nTotalwBlanks = 2*nTotal;

run_dur     = (nTotal * image_dur) + (nTotal * Blank_dur);

sz      = 768; % this is the native resolution of the MEG display
%  (actually 1024 x 768, but we restrict to a square
%  with the smaller dimension)

% contains the struct with the images from the Hermes et al 2014 paper
load('/Volumes/server/Projects/MEG/Gamma/stimuli/example_V1_electrode_faceshouses.mat');
natImages = out.image;
clear out;

scale = sz/size(natImages,1); %rescale to MEG display res
natImages = imresize(natImages,scale);

% projpth = '~/matlab/MEG_gamma';
projpth = '/Volumes/server/Projects/MEG/Gamma';

savepth = fullfile(projpth, 'stimuli/binarized');

if ~exist(savepth, 'dir'), mkdir(savepth); end

runs = 15; % Number of .mat files generated

%% pick nNatural number of natural images randomly

theseIndices = randsample(72,nNatural);


for scanNum = 1:runs
    
    %fprintf(['\n making images for scan ' num2str(scanNum) '...']);
    
    
    %% Natural Images
    
    
theseNatImages = natImages(:,:,theseIndices);

theseNatImages = repmat(theseNatImages,1,1,7);

theseNatImages = reshape(theseNatImages,768,768,7,7);

theseNatImages = permute(theseNatImages,[1 2 4 3]);

im.n = reshape(theseNatImages, 768,768,49);
    
    
%     %% white noise
%     im.wn = zeros(sz, sz, nImages, 'uint8');
%     n = 0; % white noise
%     for ii = 1:nImages
%         
%         tmp = noiseonf(sz, n);      % make a white noise stimulus       
%         im.wn(:,:,ii) = scale_images(tmp); % 8 bit integers
%         
%     end
    
%     %% binarized white noise
%     im.bin_wn = zeros(sz, sz, nImages, 'uint8');
%     n = 0; % white noise
%     for ii = 1:nImages
%         
%         tmp = noiseonf(sz, n);      % make a white noise stimulus
%         inds = tmp > median(tmp(:)); 
%         tmp(inds) = 1;
%         tmp(~inds) = 0;
%         im.bin_wn(:,:,ii) = scale_images(tmp); % 8 bit integers
%         
%     end
%     
    
    %% binarized pink noise
    im.pn = zeros(sz, sz, nImages, 'uint8');
    n = 1; % pink noise
    for ii = 1:nImages
        tmp = noiseonf(sz, n);       
        im.pn(:,:,ii) = scale_images(tmp);
        inds = tmp > median(tmp(:)); 
        tmp(inds) = 1;
        tmp(~inds) = 0;
        im.pn(:,:,ii) = scale_images(tmp); % 8 bit integers
    end
    
%     %% binarized brown noise
%     im.bn = zeros(sz, sz, nImages, 'uint8');
%     n = 2; % brown noise
%     
%     for ii = 1:nImages
%         tmp = noiseonf(sz, n);
%         tmp = tmp - min(tmp(:));
%         tmp = tmp / max(tmp(:))*255;
%         im.bn(:,:,ii) = scale_images(tmp);
%         inds = tmp > median(tmp(:)); 
%         tmp(inds) = 1;
%         tmp(~inds) = 0;
%         im.bn(:,:,ii) = scale_images(tmp); % 8 bit integers
%     end
    
    %% bars
    widths = 8;%[8 16 32 64];
    im.bars = zeros(sz,sz, length(widths) * nImages, 'uint8');
    ii = 1;
    [x, y] = meshgrid((1:sz)/sz,(1:sz)/sz);
    
    for w = 1:length(widths)
        
        for n = 1:nImages
            ph = n/nImages * 2 * pi;
            tmp = square(x*2*pi*widths(w) + ph)+1*255;
            im.bars(:,:,ii)  = scale_images(tmp);
            ii = ii+1;
        end
    end
    
%     %% Plaid
%     widths = 32;
%     im.plaid = zeros(sz,sz,  nImages, 'uint8');
%     
%     [x, y] = meshgrid((1:sz)/sz,(1:sz)/sz);
%     
%     for n = 1:nImages
%         ph = n/nImages * 2 * pi;
%         tmp = square(x*2*pi*widths + ph) + square(y*2*pi*widths(1) + ph);
%         im.plaid(:,:,n)  = scale_images(tmp);
%     end
    
    %% Grey
    im.grey = ones(sz, sz, nImages, 'uint8') * bckgrnd;
    
    %% Blanks between stimuli
    im.blank = ones(sz, sz, 1, 'uint8') * bckgrnd;
    
    
    %% concatenate all images
    % concatenates the stimfile into one 3D matrix (768x768x270)
    Images = cat(3, im.pn, im.bars, im.n, im.grey, im.blank);
    %% Creates and shuffles an index of stim conditions
    
    
    % finds the total amount of images
    nImages = size(Images, 3);
    % total amount of stimulus images
    nStim   = nImages - 1;
    % creates a randomized index
    stimindex = 1:(nStim);
    [~, randindex] = sort(rand(size(stimindex)));
    blankindex = nImages;
    
    seq = zeros(1,2*(nStim));
    seq(1:2:end) = randindex;
    seq(2:2:end) = blankindex;
    
    
  
    
    %% Categorises conditions for trigSeq
    catindex = ceil((seq)/ nImages);
    
    % catindex returns values corresponding to the category of stimuli
    % presented:

    % 1  = white noise
    % 2  = binarized white noise
    % 3  = pink noise
    % 4  = brown noise
    % 5  = gratings, width = 64
    % 6  = gratings, width = 32
    % 7  = gratings, width = 16
    % 8  = gratings, width = 8
    % 9  = plaid
    % 10 = grey
    % 11 = interspacing grey screen

    %% Photodiode Sequence
    
    diodeindex = zeros(length(catindex),1);
    for n = 1:length(catindex)
        if catindex(n) ~= max(catindex);
            diodeindex(n) = 1;
        end
    end
    
    
    
    
    % %% Preview Stimuli
    figure(2);
    for ii = 1:size(Images, 3);
        imagesc(Images(:,:,ii));
        colormap gray;
        axis image off;
        pause(0.05);
    end
    
    
    %% make stimulus parameters
    stimulus.images = Images;
    stimulus.cmap = [(1:256)' (1:256)' (1:256)'];
    stimulus.seq = seq';
    
   
    
    timing = (0:Blank_dur:run_dur);
    timing = timing(setdiff(1:length(timing), 2:3:length(timing))); 
    stimulus.seqtiming = timing(1:length(timing)-1);
    
    % Make new fixation sequence
    %   minimum time between fixation change
    min_fix_frames = round(10/image_dur); 
    %   maximum time between fixation change
    max_fix_frames = round(20/image_dur); 
    %   initialize the fixation vector with ones. some of these will change
    %   to twos.
    fix_vector = ones(size(stimulus.seq));
    counter = 1;
    
    while counter < length(fix_vector)
        % pick a random interval for this fixation where fix == 2
        this_dur = randi([min_fix_frames max_fix_frames]);
        fix_vector((1:this_dur)+counter-1) = 2;
        % pick a random interval for this fixation where fix == 1
        counter = counter + this_dur + randi([min_fix_frames max_fix_frames]);
    end
    
    % clip fixation vector in case it is longer than the image sequence
    if length(fix_vector) > length(stimulus.seq)
        fix_vector = fix_vector(1:length(stimulus.seq));
    end

    stimulus.fixSeq     = fix_vector;
    stimulus.srcRect    = [0 0 768 768];
    stimulus.destRect   = [ 128 0 896 786];
    stimulus.trigSeq    = catindex;
    stimulus.diodeSeq   = diodeindex;
    
    
    savename = fullfile(savepth, ['gammaStimuli_params' num2str(scanNum) '.mat']);
    
    save(savename, 'stimulus');
    
    fprintf('done! \n');
    
end
