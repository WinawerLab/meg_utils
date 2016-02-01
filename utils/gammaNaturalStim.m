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
% 7)whiteNoise
% 8)gratings
% 9)blank


%% parameters

background = 128; % mean pixel intensity
range = [1 225]; % pixel range
imageDuration = 1.0; % in milliseconds
blankDuration = 0.5; % ITI duration

nImages = 9; % images per category
nCategories = 9; % images excluding ITI
nTotal = nImages * nCategories;
nTotalWithITI = 2 * nTotal;

runDuration = (nTotal * imageDuration) + (nTotal * blankDuration);

size = 768; % native resolution of MEG display restricted to a square

projectPath = '/Volumes/server/Projects/MEG/Gamma';
savePath = fullfile(projectPath, 'stimuli/natural_images');
imagePath = fullfile(projectPath, '/natural_images/nat_images_before');
if ~exist(savePath, 'dir'), mkdir(savePath); end

% number of .mat files generated
runs = 15;

%% houses

houseIndices = [1 14 72]; % file number of selected image

house1 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(1))));
house2 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(2))));
house3 = imread(fullfile(imagePath, sprintf('/nat_image%d.png', houseIndices(3))));

houseH = zeros(size, size, nImages);
houseM = zeros(size, size, nImages);
houseL = zeros(size, size, nImages);

% high contrast
for i = 1:3
    houseH(:,:,i) = normalizeLaplacian(house1, 'high');
    houseH(:,:,i+3) = normalizeLaplacian(house2, 'high');
    houseH(:,:,i+6) = normalizeLaplacian(house3, 'high');
end

% medium contrast
for i = 1:3
    houseH(:,:,i) = normalizeLaplacian(house1, 'high');
    houseH(:,:,i+3) = normalizeLaplacian(house2, 'high');
    houseH(:,:,i+6) = normalizeLaplacian(house3, 'high');
end


for i = 1:9
    imshow(houseH(:,:,i));
    pause(0.5);

end
