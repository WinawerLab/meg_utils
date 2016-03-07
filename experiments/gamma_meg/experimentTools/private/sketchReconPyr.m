function reconIm = sketchReconPyr(pyramid, LPR)
% SKETCH RECON PYR - a teaching example of how to
% reconstruct an image from a Laplacian pyramid, Dec 2015

layers = length(pyramid);

lowPass = LPR;
for ii = layers:-1:1
    expandLP = imresize(lowPass, 2);
    lowPass = expandLP + pyramid{ii};
end
reconIm = lowPass;

end

