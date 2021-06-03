function [Masks] = segmentThreshold(I)

% adaptive thresholding
T = adaptthresh(I);
Masks = imbinarize(I, T);

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
Masks = imclearborder(Masks);

end