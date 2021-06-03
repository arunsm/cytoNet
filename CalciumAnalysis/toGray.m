function im = toGray(im)
if size(im, 3) >= 3
    im = rgb2gray(im(:, :, 1:3));
else
    im = im(:, :, 1);
end
end
    