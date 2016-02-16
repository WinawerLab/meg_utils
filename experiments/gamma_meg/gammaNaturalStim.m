%% gammaNaturalStimuli
%
% generates and saves stimulus strucutres for running meg natural image
% experiment with vistadisp. 
%
% The files that are generated contain the following fields:
%   images (X x Y X n) array (X = rows, Y = cols, n = number of images)
%   cmap: 256 x 3 colormap (each column is 1:256 unless we are doing
%               something unusual with images)
%   seq: numframes x 1 indicating the sequence of images to be presented
%               during the experiment
%   seqtiming: numframes x 1 indicating the time (in seconds) at which each image is shown
%   fixSeq:    numframes x 1 indicating the fixation image to be shown
%                   (typically 1s and 2s that slowly alternate)
%
%
%  Stimulus set contains faces, houses, gratings, and noise at 3 different
%  contrast levels
% 
% stimulus categories:
% 1)houseH
% 2)houseM
% 3)houseL
% 4)faceH
% 5)faceM
% 6)faceL
% 7)pinkNoise
% 8)gratings
% 9)blank
%
% Nicholas Chua (2015), script for generating noise and gratings by Dora Hermes


%% parameters
saveFiles = false;
visualizeImages = true;


background = 128; % mean pixel intensity
range = [1 255]; % pixel range
imageDuration = 1.0; % in milliseconds
blankDuration = 0.5; % ITI duration

% Each run will have 9 repeats of 13 images = 117 trials
nImages         = 9;    % images per category
nCategories     = 13;   % images excluding ITI
nTotal          = nImages * nCategories;
nTotalWithITI   = 2 * nTotal;

runDuration = (nTotal * imageDuration) + (nTotal * blankDuration);

sz = 768; % native resolution of MEG display restricted to a square

% make a circular mask
[x, y] = meshgrid(linspace(-1,1,sz));
mask = x.^2 + y.^2 <= .8;
mask = true(sz);

projectPath = '/Volumes/server/Projects/MEG/Gamma';
savePath    = fullfile(projectPath, 'stimuli/natural_images');
imagePath   = fullfile(projectPath, 'cerebral_cortex_datashare/example_V1_electrode_faceshouses.mat');
if ~exist(savePath, 'dir'), mkdir(savePath); end

addpath(genpath('~/matlab/git/vistadisp/'));

scale_images = @(x) uint8((x - min(x(:))) / (max(x(:)) - min(x(:))) * diff(range) + min(range));

% number of .mat files generated - each one has a different random order of
% trials
totalRuns = 1;



%% houses
nHouseImages = 3; % number of different house images
% File number of selected image. 
%   Natural scene images come from Hermes et al, 2014. These numbers are
%   indices into the order of these images
houseIndices = [1 14 72]; 

load(imagePath, 'out');
imagesOrig = out.image+128; clear out;
house1 = imresize( imagesOrig(:,:, houseIndices(1)), [sz sz]);
house2 = imresize( imagesOrig(:,:, houseIndices(2)), [sz sz]);
house3 = imresize( imagesOrig(:,:, houseIndices(3)), [sz sz]);

houseH = zeros(sz, sz, nImages/nHouseImages);
houseM = zeros(sz, sz, nImages/nHouseImages);
houseL = zeros(sz, sz, nImages/nHouseImages);

% high contrast
houseH(:,:,1) = normalizeLaplacian(house1, 'high', mask);
houseH(:,:,2) = normalizeLaplacian(house2, 'high', mask);
houseH(:,:,3) = normalizeLaplacian(house3, 'high', mask);




% medium contrast
houseM(:,:,1) = normalizeLaplacian(house1, 'medium', mask);
houseM(:,:,2) = normalizeLaplacian(house2, 'medium', mask);
houseM(:,:,3) = normalizeLaplacian(house3, 'medium', mask);

% low contrast
houseL(:,:,1) = normalizeLaplacian(house1, 'low', mask);
houseL(:,:,2) = normalizeLaplacian(house2, 'low', mask);
houseL(:,:,3) = normalizeLaplacian(house3, 'low', mask);


% repeat images to obtain nImages number of houses
houseH = repmat(houseH, 1, 1, nImages/nHouseImages);
houseM = repmat(houseM, 1, 1, nImages/nHouseImages);
houseL = repmat(houseL, 1, 1, nImages/nHouseImages);

% scale and convert to uint8
% for i = 1:size(houseH, 3)
%     houseH(:,:,i) = scale_images(houseH(:,:,i));
%     houseM(:,:,i) = scale_images(houseM(:,:,i));
%     houseL(:,:,i) = scale_images(houseL(:,:,i));
% end
houseH = uint8(houseH);
houseM = uint8(houseM);
houseL = uint8(houseL);


stdH = std(double(reshape(houseH, sz*sz,[])));
stdM = std(double(reshape(houseM, sz*sz,[])));
stdL = std(double(reshape(houseL, sz*sz,[])));

figure, 
for ii = 1:3
    subplot(3,3,0+ii); histogram(houseH(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdH(ii)));
    
    subplot(3,3,3+ii); histogram(houseM(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdM(ii)));
    
    subplot(3,3,6+ii); histogram(houseL(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdL(ii)));
end

%% faces
nFaceImages = 3; % number of different house images
faceIndices = [16 39 66]; % file number of selected image

face1 = imresize( imagesOrig(:,:, faceIndices(1)), [sz sz]);
face2 = imresize( imagesOrig(:,:, faceIndices(2)), [sz sz]);
face3 = imresize( imagesOrig(:,:, faceIndices(3)), [sz sz]);

% faceH = uint8(zeros(sz, sz, nImages/nFaceImages));
% faceM = uint8(zeros(sz, sz, nImages/nFaceImages));
% faceL = uint8(zeros(sz, sz, nImages/nFaceImages));

faceH = zeros(sz, sz, nImages/nFaceImages);
faceM = zeros(sz, sz, nImages/nFaceImages);
faceL = zeros(sz, sz, nImages/nFaceImages);

% high contrast 
faceH(:,:,1) = normalizeLaplacian(face1, 'high', mask);
faceH(:,:,2) = normalizeLaplacian(face2, 'high', mask);
faceH(:,:,3) = normalizeLaplacian(face3, 'high', mask);

% medium contrast
faceM(:,:,1) = normalizeLaplacian(face1, 'medium', mask);
faceM(:,:,2) = normalizeLaplacian(face2, 'medium', mask);
faceM(:,:,3) = normalizeLaplacian(face3, 'medium', mask);

% low contrast
faceL(:,:,1) = normalizeLaplacian(face1, 'low', mask);
faceL(:,:,2) = normalizeLaplacian(face2, 'low', mask);
faceL(:,:,3) = normalizeLaplacian(face3, 'low', mask);

% repeat images to obtain nImages number of face
faceH = repmat(faceH, 1, 1, nImages/nFaceImages);
faceM = repmat(faceM, 1, 1, nImages/nFaceImages);
faceL = repmat(faceL, 1, 1, nImages/nFaceImages);

% scale and convert to uint8
% for i = 1:size(faceH, 3)
%     faceH(:,:,i) = scale_images(faceH(:,:,i));
%     faceM(:,:,i) = scale_images(faceM(:,:,i));
%     faceL(:,:,i) = scale_images(faceL(:,:,i));
% end

faceH = uint8(faceH); 
faceM = uint8(faceM);
faceL = uint8(faceL);


stdH = std(double(reshape(faceH, sz*sz,[])));
stdM = std(double(reshape(faceM, sz*sz,[])));
stdL = std(double(reshape(faceL, sz*sz,[])));

figure, 
for ii = 1:3
    subplot(3,3,0+ii); histogram(faceH(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdH(ii)));
    
    subplot(3,3,3+ii); histogram(faceM(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdM(ii)));
    
    subplot(3,3,6+ii); histogram(faceL(:,:,ii)); xlim([0 255]);
    title(sprintf('Std = %6.2f', stdL(ii)));
end



%% 
for run = 1:totalRuns
    fprintf(['\n making images for scan ' num2str(run) '...']);
    
    %% pink noise
    
    pinkNoise  = zeros(sz, sz, nImages*3, 'uint8');
    pinkNoiseH = zeros(sz, sz, nImages, 'uint8');
    pinkNoiseM = zeros(sz, sz, nImages, 'uint8');
    pinkNoiseL = zeros(sz, sz, nImages, 'uint8');
    
    n = 1; % pink noise
    
    for ii = 1:nImages*3 % 3 contrast conditions
        tmp = noiseonf(sz, n);
        inds = tmp > median(tmp(:));
        tmps(inds) = 1;
        tmps(~inds) = 0;
        pinkNoise(:,:,ii) = scale_images(tmp);
    end
    
    for iii = 1:nImages
        pinkNoiseH(:,:,iii) = normalizeLaplacian(pinkNoise(:,:,iii),'high');
        pinkNoiseM(:,:,iii) = normalizeLaplacian(pinkNoise(:,:,iii+9),'medium');
        pinkNoiseL(:,:,iii) = normalizeLaplacian(pinkNoise(:,:,iii+18),'low');
    end
    
    
    
    %% gratings
    
    cycles = 64; % [8 16 32 64]; cycles per image
    gratings = zeros(sz, sz, nImages * length(cycles) * 3);

    [x, y] = meshgrid((1:sz)/sz, (1:sz)/sz);
    
    for iii = 1:length(cycles)
        for n = 1:nImages*3 % 3 contrast conditions
            ph = n/nImages * 2 * pi;
            tmp = square(x * 2 * pi * cycles(iii) + ph) + 1 * 255;
            gratings(:,:,n) = scale_images(tmp);
        end
    end
    
    for iii = 1:nImages
        gratingsH(:,:,iii) = normalizeLaplacian(gratings(:,:,iii),'high');
        gratingsM(:,:,iii) = normalizeLaplacian(gratings(:,:,iii+9),'medium');
        gratingsL(:,:,iii) = normalizeLaplacian(gratings(:,:,iii+18),'low');
    end
    
    %% blank stimuli
    
    blank = ones(sz, sz, nImages, 'uint8') * background;
    
    %% intertrail interval
    
    ITI = ones(sz, sz, 1, 'uint8') * background;
    
    %% Concatenate all stimuli images and shuffle
    
    images = cat(3, houseH, houseM, houseL, faceH, faceM, faceL,...
        pinkNoiseH, pinkNoiseM, pinkNoiseL, gratingsH, gratingsM,...
        gratingsL, blank, ITI);
    
    stimIndex = 1:(nTotal); % index of all stimuli images
    [~, randIndex] = sort(rand(size(stimIndex))); % shuffle index randomly
    ITIindex = size(images, 3); % ITI is the last image in the concatenated array
    
    sequence = zeros(1, 2 * nTotal); % every other image is an ITI
    sequence(1:2:end) = randIndex;
    sequence(2:2:end) = ITIindex;
    
    % index of stimulus categories
    categoryIndex = ceil((sequence/nImages));
    
    stdev = NaN(1,size(images,3));
    mn    = NaN(1,size(images,3));
    % apply soft circular aperture
    for i = 1:size(images,3)
        [images(:,:,i), stdev(i), mn(i)] = cosineMask(images(:,:,i));
    end
    
    %% save images
%     sampleSavePath = fullfile(savePath, 'sampleImages');
%     
%     for i = 1:nImages:size(images,3)
%         saveName = fullfile(sampleSavePath, ['stimImage' num2str(i) '.png']);
%         imwrite(images(:,:,i), saveName);
%     end
%     
%     for i = 2:nImages:size(images,3)
%         saveName = fullfile(sampleSavePath, ['stimImage' num2str(i) '.png']);
%         imwrite(images(:,:,i), saveName);
%     end
    %% photodiode
    
    diodeIndex = zeros(length(categoryIndex), 1);
    for i = 1:length(categoryIndex)
        if categoryIndex(i) ~= max(categoryIndex)
            diodeIndex(i) = 1;
        end
    end
    
    %% preview stimuli
    if visualizeImages
        figure;
        for ii = 1:size(images, 3);
            %imagesc(images(:,:,ii),range);
            %colormap gray;
            %axis image off;
            imshow(images(:,:,ii));
            title(sprintf('Mean = %6.2f; Std = %6.2f', mn(ii), stdev(ii)));
            waitforbuttonpress;%  pause(0.2);
        end
    end
    %% stimulus parameters
    
    stimulus.images = images;
    stimulus.cmap = [(1:256)' (1:256)' (1:256)'];
    stimulus.seq = sequence';
    
    % timing
    timing = (0:blankDuration:runDuration);
    timing = timing(setdiff(1:length(timing), 2:3:length(timing)));
    stimulus.seqtiming = timing(1:length(timing) - 1);
    
    % fixation sequence
    minFixFrames = round(10/imageDuration);
    maxFixFrames = round(20/imageDuration);
    fixVector = ones(size(stimulus.seq));
    counter = 1;
    
    while counter < length(fixVector)
        this_dur = randi([minFixFrames maxFixFrames]); % random duration
        fixVector((1:this_dur)+counter -1) = 2;
        counter = counter + this_dur + randi([minFixFrames maxFixFrames]);
    end
    
    if length(fixVector) > length(stimulus.seq) % clip if longer than run
        fixVector = fixVector(1:length(stimulus.seq));
    end
    
    stimulus.fixSeq = fixVector;
    stimulus.srcRect = [0 0 768 768];
    stimulus.destRect = [128 0 896 786];
    stimulus.trigSeq = categoryIndex;
    stimulus.diodeSeq = diodeIndex;
    
    if saveFiles
        saveName = fullfile(savePath, ['gammaStimuli_params' num2str(run) '.mat']);
        
        save(saveName, 'stimulus');
    end
end




