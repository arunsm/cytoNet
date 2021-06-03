%This function takes an image file and segments it though several steps.
%First, the object is smoothed, and a seed mask is created from the
%regional maximums. Then, the image is segmented, and the segmentation is
%split by a watershed transform of the seed mask. 

%Input: name of image to be segmented
function [maxImage, segmented] = segmentNPC(filename)

info = imfinfo(filename);
mask = zeros(info(1).Height, info(1).Width);
imgs = zeros(info(1).Height, info(1).Width, numel(info));
maxImage = zeros(size(mask));
for i = 1:numel(info)
    im = im2double(toGray(imread(filename,i))); %Converts image to double type for uniform computation
    imgs(:, :, i) = im;
    maxImage = max(im, maxImage);

    %Perform initial segmentation using Otsu's Method.
    thresh = graythresh(im);
    bw = im2bw(im, thresh);
    
    %Clean up stray pixels and neurites using opening with a disk of radius 2.
    bw = imopen(bw, strel('disk', 5));
    
    
    %Add newest mask to overall mask
    mask = mask | bw;
    %fprintf('%i/%i\n', i, numel(info));
end

%Use watershed to split connected cells
%im = toGray(imread(filename));
im = mean(imgs, 3); 
im = mat2gray(im);
im2 = imcomplement(im); %New pixel = 1 - old pixel
im3 = imhmin(im2, .08); %Supresses minima less less than .2 to prevent oversegmentation
L = watershed(im3); %Computes "catchbasin" regions; goal is to find "ridgelines" where catchbasins meet (presumably where objects in image are touching)
mask(L == 0) = 0; %Ridgelines are on L == 0

%Return images
mask = imfill(mask, 'holes');
segmented = mask;
maxImage = mat2gray(maxImage);
end

