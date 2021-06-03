
function returnStruct = readMasks(fileName)

returnStruct = struct;

imageInfo = imfinfo(fileName);

imageWidth = imageInfo.Width; % image width
imageHeight = imageInfo.Height; % image height
nLayers = size(imageInfo, 1);

cc_cb = struct;

cc_cb.Connectivity = 8; % set randomly - not relevant in this case
cc_cb.ImageSize = [imageHeight imageWidth];
cc_cb.NumObjects = nLayers;
cc_cb.PixelIdxList = {};

allMasks = false(imageHeight, imageWidth);

masks = cell(nLayers, 1);
MaskBoundaries = cell(nLayers, 1);

for i = 1:nLayers
    mask_i = imread(fileName, i); % reading i'th layer of tif image
    masks{i} = mask_i;
    
    MaskBoundariesCell = bwboundaries(mask_i, 'noholes');
    MaskBoundaries{i} = MaskBoundariesCell{1};
    
    allMasks = allMasks | mask_i;
    
    [pixel_idx_i] = find(mask_i); % linear pixel indices of i'th mask
    cc_cb.PixelIdxList{i} = pixel_idx_i;
end

returnStruct.MasksCellArray = masks;
returnStruct.CellBodyStruct = cc_cb;
returnStruct.MaskBoundaries = MaskBoundaries;

end