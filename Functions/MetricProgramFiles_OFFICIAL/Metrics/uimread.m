function I = uimread(filename, index)
finfo = imfinfo(filename);
% Some image file formats do not support multiple images per file.  For the
% first and perhaps only image in a file, call imread without an index
% argument.
if index == 1
    [I map] = imread(filename);
else
    [I map] = imread(filename, index);
end

% Convert color and indexed images to grayscale
colortp = finfo(index).ColorType;
switch colortp
    case 'truecolor'
        I = rgb2gray(I);
    case 'indexed'
        I = ind2rgb(I, map);
        I = rgb2gray(I);
end

% Convert grayscale values to type double ranging between 0 and 1
I = im2double(I);

%switch class(I)
%    case 'single'
%        max = 1;
%    case 'double'
%        max = 1;
%    case 'uint8'
%        max = (2^8) - 1;
%    case 'uint16'
%        max = (2^16) - 1;
%    case 'uint32'
%        max = (2^32) - 1;
%    case 'uint64'
%        max = (2^64) - 1;
%    otherwise
%       error('uimread:unexpected', 'Unexpected image class: %s', class(I));
%end
%I = mat2gray(I, [0 max]);
end

