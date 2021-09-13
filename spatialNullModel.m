function [spatialNullMaskDilated, spatialRandomAdjacencyMatrix, spatialRandomCellLocations] = spatialNullModel(mask, threshold, adjacencyType)

%% create random spatial map of cells
CC = bwconncomp(mask);
nNodes = CC.NumObjects;

% compute average equivalent diameter of cells
d = regionprops(CC, 'EquivDiameter');
diameters = cat(1, d.EquivDiameter);
avgeDiameter = mean(diameters);

imageDimensions = CC.ImageSize;

% place dots at random X and Y locations
randLocationsX = randperm(imageDimensions(1), nNodes)';
randLocationsY = randperm(imageDimensions(2), nNodes)';
spatialRandomCellLocations = [randLocationsY, randLocationsX];

spatialNullMask = false(size(mask));
for i = 1:nNodes
    spatialNullMask(randLocationsX(i), randLocationsY(i)) = true;
end

% dilate mask with sphere of diameter equal to the equivalent diameter
% of cells
spatialNullMaskDilated = imdilate(spatialNullMask, strel('sphere', floor(avgeDiameter/2)));

%% compute adjacency matrix
if adjacencyType == 1
    % Splitting each cell mask into a separate image in a cell array
    masksCellArray = cell(nNodes, 1);
    for i = 1:nNodes
        mask_i = false(imageDimensions);
        mask_i(CC.PixelIdxList{i}) = true;
        masksCellArray{i} = mask_i;
    end
    
    allMasksExtended = cell([nNodes 1]); % initialize new binary image containing all extended masks
    filter = [0 1 0; 1 1 1; 0 1 0]; % filter for mask extension
    
    % extending masks out by 2 pixels
    for i = 1:nNodes
        mask_i = masksCellArray{i};
        mask_i_extended_1pixel = logical(imfilter(mask_i, filter)); % extending mask i based on filter
        mask_i_extended_2pixels = logical(imfilter(mask_i_extended_1pixel, filter));
        
        allMasksExtended{i} = mask_i_extended_2pixels;
    end
    
    borderLength = squareform(sharedPixels(allMasksExtended));
    spatialRandomAdjacencyMatrix = double(borderLength > threshold);
end

if adjacencyType == 2
    diameters = repmat(avgeDiameter, nNodes, 1);
    distanceThreshold = threshold*squareform(pdist(diameters, @averageDiameter));
    weights = squareform(pdist(spatialRandomCellLocations, 'Euclidean'));
    spatialRandomAdjacencyMatrix = double(weights < distanceThreshold);
end

%% optional display of original graph and spatial random graph
%f = figure('Visible', 'Off');
% figure; ax = gca;
% imshow(spatialNullMaskDilated, 'Border', 'Tight'); hold(ax, 'on');
% 
% if ~isempty(spatialRandomAdjacencyMatrix)
%     gplot(spatialRandomAdjacencyMatrix, spatialRandomCellLocations, 'y');
% end
end