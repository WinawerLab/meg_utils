function [pyramid, LPR] = sketchBuildPyr(im, layers)
% SKETCH BUILD PYRAMID - a teaching example of how to
% build a Laplacian pyramid, Dec 2015
    pyramid = cell(1, layers);

    lowPass = im;
    for ii = 1:layers
        reduced = imresize(lowPass, 0.5);
        expanded = imresize(reduced, 2);
        diff = lowPass - expanded;

        lowPass = reduced;
        pyramid{ii} = diff;
    end
    LPR = lowPass;
end

