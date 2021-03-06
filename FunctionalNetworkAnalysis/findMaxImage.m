% Function to find maximum intensity projection of a stack of images

function maxImage = findMaxImage(filename)

info = imfinfo(filename);
maxImage = zeros(info(1).Height, info(1).Width);
for i = 1:numel(info)
    im = im2double(imread(filename,i)); % convert image to double type for uniform computation
    maxImage = max(im, maxImage);
end

end