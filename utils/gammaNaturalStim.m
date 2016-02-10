%% gammaNaturalStimuli
% generates .mat files containing the stimuli images
% contains faces and houses of different contrast levels
% gratings and high contrast noise
% 
% stimulus categories:
% 1)houseH
% 2)houseM
% 3)houseL
% 4)faceH
% 5)faceM
% 6)faceL
% 7)binarizedWhiteNoise
% 8)gratings
% 9)blank
%
% Nicholas Chua (2015), script for generating noise and gratings by Dora Hermes


%% parameters

background = 128; % mean pixel intensity
range = [1 225]; % pixel range
imageDuration = 1.0; % in milliseconds
blankDuration = 0.5; % ITI duration

nImages = 9; % images per category
nCategories = 13; % images excluding ITI
nTotal = nImages * nCategories;
nTotalWithITI = 2 * nTotal;

runDuration = (nTotal * imageDuration) + (nTotal * blankDuration);

sz = 768; % native resolution of MEG display restricted to a square

projectPath = '/Volumes/server/Projects/MEG/Gamma';
savePath = fullfile(projectPath, 'stimuli/natural_images');
imagePath = fullfile(projectPath, '/natural_images_tools/nat_images_before');
if ~exist(savePath, 'dir'), mkdir(savePath); end

addpath(genpath('~/matlab/git/meg_utils'));

addpath(genpath('~/matlab/git/vistadisp/'));

scale_images = @(x) uint8((x - min(x(:))) / (max(x(:)) - min(x(:))) * diff(range) + min(range));

% number of .mat files generated
totalRuns = 12;

saveFiles = true;
visualizeImages = true;

%% houses
nHouseImages = 3; % number of different house images
houseIndices = [1 14 72]; % file number of selected image

house1 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(1))));
house2 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(2))));
house3 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(3))));

% houseH = uint8(zeros(sz, sz, nImages/nHouseImages));
% houseM = uint8(zeros(sz, sz, nImages/nHouseImages));
% houseL = uint8(zeros(sz, sz, nImages/nHouseImages));

houseH = zeros(sz, sz, nImages/nHouseImages);
houseM = zeros(sz, sz, nImages/nHouseImages);
houseL = zeros(sz, sz, nImages/nHouseImages);

% high contrast
houseH(:,:,1) = normalizeLaplacian(house1, 'high');
houseH(:,:,2) = normalizeLaplacian(house2, 'high');
houseH(:,:,3) = normalizeLaplacian(house3, 'high');

% medium contrast
houseM(:,:,1) = normalizeLaplacian(house1, 'medium');
houseM(:,:,2) = normalizeLaplacian(house2, 'medium');
houseM(:,:,3) = normalizeLaplacian(house3, 'medium');

% low contrast
houseL(:,:,1) = normalizeLaplacian(house1, 'low');
houseL(:,:,2) = normalizeLaplacian(house2, 'low');
houseL(:,:,3) = normalizeLaplacian(house3, 'low');

% repeat images to obtain nImages number of houses
houseH = repmat(houseH, 1, 1, nImages/nHouseImages);
houseM = repmat(houseM, 1, 1, nImages/nHouseImages);
houseL = repmat(houseL, 1, 1, nImages/nHouseImages);

% scale and convert to uint8
for i = 1:size(houseH, 3)
    houseH(:,:,i) = scale_images(houseH(:,:,i));
    houseM(:,:,i) = scale_images(houseM(:,:,i));
    houseL(:,:,i) = scale_images(houseL(:,:,i));
end

%% faces
nFaceImages = 3; % number of different house images
faceIndices = [16 39 66]; % file number of selected image

face1 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', faceIndices(1))));
face2 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', faceIndices(2))));
face3 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', faceIndices(3))));

% faceH = uint8(zeros(sz, sz, nImages/nFaceImages));
% faceM = uint8(zeros(sz, sz, nImages/nFaceImages));
% faceL = uint8(zeros(sz, sz, nImages/nFaceImages));

faceH = zeros(sz, sz, nImages/nFaceImages);
faceM = zeros(sz, sz, nImages/nFaceImages);
faceL = zeros(sz, sz, nImages/nFaceImages);

% high contrast
faceH(:,:,1) = normalizeLaplacian(face1, 'high');
faceH(:,:,2) = normalizeLaplacian(face2, 'high');
faceH(:,:,3) = normalizeLaplacian(face3, 'high');

% medium contrast
faceM(:,:,1) = normalizeLaplacian(face1, 'medium');
faceM(:,:,2) = normalizeLaplacian(face2, 'medium');
faceM(:,:,3) = normalizeLaplacian(face3, 'medium');

% low contrast
faceL(:,:,1) = normalizeLaplacian(face1, 'low');
faceL(:,:,2) = normalizeLaplacian(face2, 'low');
faceL(:,:,3) = normalizeLaplacian(face3, 'low');

% repeat images to obtain nImages number of face
faceH = repmat(faceH, 1, 1, nImages/nFaceImages);
faceM = repmat(faceM, 1, 1, nImages/nFaceImages);
faceL = repmat(faceL, 1, 1, nImages/nFaceImages);

% scale and convert to uint8
for i = 1:size(faceH, 3)
    faceH(:,:,i) = scale_images(faceH(:,:,i));
    faceM(:,:,i) = scale_images(faceM(:,:,i));
    faceL(:,:,i) = scale_images(faceL(:,:,i));
end

%% 
for run = 1:totalRuns
    fprintf(['\n making images for scan ' num2str(run) '...']);
    
    %% white noise
    
    whiteNoise = zeros(sz, sz, nImages*3, 'uint8');
    whiteNoiseH = zeros(sz, sz, nImages, 'uint8');
    whiteNoiseM = zeros(sz, sz, nImages, 'uint8');
    whiteNoiseL = zeros(sz, sz, nImages, 'uint8');
    
    n = 0; % white noise
    
    for ii = 1:nImages*3 % 3 contrast conditions
        tmp = noiseonf(sz, n);
        inds = tmp > median(tmp(:));
        tmps(inds) = 1;
        tmps(~inds) = 0;
        whiteNoise(:,:,ii) = scale_images(tmp);
    end
    
    for iii = 1:nImages
        whiteNoiseH(:,:,iii) = normalizeLaplacian(whiteNoise(:,:,iii),'high');
        whiteNoiseM(:,:,iii) = normalizeLaplacian(whiteNoise(:,:,iii+9),'medium');
        whiteNoiseL(:,:,iii) = normalizeLaplacian(whiteNoise(:,:,iii+18),'low');
    end
    
    
    
    %% gratings
    
    widths = 8; % [8 16 32 64]; pixles
    gratings = zeros(sz, sz, nImages * length(widths) * 3);
    [x, y] = meshgrid((1:sz)/sz, (1:sz)/sz);
    
    for iii = 1:length(widths)
        for n = 1:nImages*3 % 3 contrast conditions
            ph = n/nImages * 2 * pi;
            tmp = square(x * 2 * pi * widths(iii) + ph) + 1 * 255;
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
        whiteNoiseH, whiteNoiseM, whiteNoiseL, gratingsH, gratingsM,...
        gratingsL, blank, ITI);
    
    stimIndex = 1:(nTotal); % index of all stimuli images
    [~, randIndex] = sort(rand(size(stimIndex))); % shuffle index randomly
    ITIindex = size(images, 3); % ITI is the last image in the concatenated array
    
    sequence = zeros(1, 2 * nTotal); % every other image is an ITI
    sequence(1:2:end) = randIndex;
    sequence(2:2:end) = ITIindex;
    
    % index of stimulus categories
    categoryIndex = ceil((sequence/nImages));
    
    
    % apply soft circular aperture
    for i = 1:size(images,3)
        images(:,:,i) = cosineMask(images(:,:,i));
    end
    
    %% save images
    sampleSavePath = fullfile(savePath, 'sampleImages');
    
    for i = 1:nImages:size(images,3)
        saveName = fullfile(sampleSavePath, ['stimImage' num2str(i) '.png']);
        imwrite(images(:,:,i), saveName);
    end
    
    for i = 2:nImages:size(images,3)
        saveName = fullfile(sampleSavePath, ['stimImage' num2str(i) '.png']);
        imwrite(images(:,:,i), saveName);
    end
    %% photodiode
    
    diodeIndex = zeros(length(categoryIndex), 1);
    for i = 1:length(categoryIndex)
        if categoryIndex(i) ~= max(categoryIndex)
            diodeIndex(i) = 1;
        end
    end
    
    %% preview stimuli
    if visualizeImages
        figure(2);
        for ii = 1:size(images, 3);
            imagesc(images(:,:,ii),range);
            colormap gray;
            axis image off;
            pause(0.2);
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




