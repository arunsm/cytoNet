% Function to compute a spatial adjacency matrix based on border pixel
% overlap (type I graph) or centroid-centroid proximity (type II graph)
% from a binary mask
%
% INPUT: 
% mask - binary mask defining ROIs
% threshold - scalar representing border pixel overlap threshold for type I
% graphs and scaling factor for type II graphs
% adjacencyType - indicator for type of graph, adjacencyType = 1 will
% create a border pixel overlap graph (type I) and adjacencyType = 2 will
% create a centroid-centroid proximity graph (type II)
%
% OUTPUT: cellInfoAllCells - MATLAB structure containing graph information
% for current image

function cellInfoAllCells = generateSpatialGraph(mask, threshold, adjacencyType)

% creating a connected component structure to represent masks
CC = bwconncomp(mask);
nNodes = CC.NumObjects;
imageDimensions = CC.ImageSize;

% Extracting cell centroids
c = regionprops(CC, 'Centroid');
cellLocations = cat(1, c.Centroid);

% Computing equivalent diameters of cells
d = regionprops(CC, 'EquivDiameter');
diameters = cat(1, d.EquivDiameter);

cellInfoAllCells = struct;
cellInfoAllCells.mask = mask;
cellInfoAllCells.nNodes = nNodes;
cellInfoAllCells.cellLocations = cellLocations;
cellInfoAllCells.graphTypeTag = 'spatial';
cellInfoAllCells.threshold = threshold;
cellInfoAllCells.adjacencyType = adjacencyType;

%% Adjacency Matrix based on shared border pixels

if adjacencyType == 1
    % Splitting each cell mask into a separate image in a cell array
    MasksCellArray = cell(nNodes, 1);
    for i = 1:nNodes
        mask_i = false(imageDimensions);
        mask_i(CC.PixelIdxList{i}) = true;
        MasksCellArray{i} = mask_i;
    end
    
    allMasksExtended = cell([nNodes 1]); % initialize new binary image containing all extended masks
    filter = [0 1 0; 1 1 1; 0 1 0]; % filter for mask extension
    
    % extending masks out by 2 pixels
    for i = 1:nNodes
        mask_i = MasksCellArray{i};
        mask_i_extended_1pixel = logical(imfilter(mask_i, filter)); % extending mask i based on filter
        mask_i_extended_2pixels = logical(imfilter(mask_i_extended_1pixel, filter));
        
        allMasksExtended{i} = mask_i_extended_2pixels;
    end
    
    borderLength = squareform(sharedPixels(allMasksExtended));
    adjacencyMatrixBinary = double(borderLength > threshold);
end

%% Adjacency Matrix based on centroid-centroid proximity

if adjacencyType == 2
    % Setting distance threshold for each pair of nodes as the average
    % diameter of that pair of cells multiplied by scaling factor
    % (threshold)
    distanceThreshold = threshold*squareform(pdist(diameters, @averageDiameter));
    
    %distanceThreshold = 282; % hard threshold option?
    
    % Computing adjacency matrix
    weights = squareform(pdist(cellLocations, 'Euclidean'));
    adjacencyMatrixBinary = double(weights < distanceThreshold);
end

% set diagonal elements of Adjacency Matrix to zero
adjacencyMatrixBinary = adjacencyMatrixBinary - diag(diag(adjacencyMatrixBinary));

cellInfoAllCells.adjacencyMatrixBinary = adjacencyMatrixBinary;
end