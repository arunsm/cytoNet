function [Masks] = segmentWatershed(I)

% adaptive thresholding
T = adaptthresh(I);
Masks = imbinarize(I, T);

Masks = imclearborder(Masks);

% performing median filtering to eliminate salt and pepper noise
Masks = medfilt2(Masks, [5 5]);

% using 20th percentile of object sizes to filter out noise
% CC = bwconncomp(Mask);
% Prop = regionprops(CC, 'Area');
% ObjectSizes = cat(1, Prop.Area);
% ObjectFilterSize = prctile(ObjectSizes, 20);
% Mask = bwareaopen(Mask, ObjectFilterSize);

% eliminating objects under 50 pixels
Masks = bwareaopen(Masks, 50);
Masks = imfill(Masks, 'holes');

% watershed seed generation using regional minima in distance transform
D = -bwdist(~Masks);
LocalMinima = imextendedmin(D, 1);
D2 = imimposemin(D, LocalMinima);
L = watershed(D2);
Masks(L == 0) = 0;

Masks = imclearborder(Masks);
end
