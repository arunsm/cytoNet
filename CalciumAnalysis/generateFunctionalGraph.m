% Function to extract calcium time series data from an image stack and
% compute a functional graph based on correlations among time series
%
% Input: filename - stacked tif image file you wish to generate a graph for
%        timetotal - total timeframe image set covers
%        timestop - timepoint you wish to stop reading in data
%        mask (optional) - mask defining ROIs for time series extraction;
%        if no mask is provided, the default segmentation program is run on
%        the maximum intensity projection image
% Output: cellInfo - cell array containing time series data for all active
%        cells (those displaying >3 transients)
%        cellInfoAllCells - cell array containing time series data for all 
%        cells
%        mask - mask used for computation
%        maxImage - maximum intensity projection image
%        functionalAdjacencyMatrix - adjacency matrix specifying
%        connections among ROIs based on time series correlations

function [errorReport, cellInfo, cellInfoAllCells, mask, maxImage, functionalAdjacencyMatrix] = generateFunctionalGraph(filename, timetotal, timestop, varargin)

%% Setting up
info = imfinfo(filename);
maxImage = findMaxImage(filename); % extract maximum intensity projection image

% if mask not provided as input, run default segmentation - NEED TO FIX
if nargin < 4
    mask = segmentWatershed(maxImage);
end

framestop = round(numel(info)*(timestop/timetotal));
mask = imresize(mask, [info(1).Height, info(1).Width]); % ensure mask and images are same dimensions
CC = bwconncomp(mask);
nCells = CC.NumObjects;

if nCells == 0
    fn = getFileName(filename);
    errorReport = makeErrorStruct(['no cells found in file ', fn], 1);
    cellInfo = [];
    cellInfoAllCells = [];
    maxImage = [];
    functionalAdjacencyMatrix = [];
    return;
end

%% Gathering info on ROIs
nActiveCells = CC.NumObjects;
cellInfoAllCells = zeros(nActiveCells, 6, framestop);
stats = regionprops(CC, 'Centroid');

%% gathering (x, y) coordinates and cell index
for k = 1:CC.NumObjects
    current = stats(k);
    location = current.Centroid;
    cellInfoAllCells(k, 1, :) = location(1); % storing x location
    cellInfoAllCells(k, 2, :) = location(2); % storing y location
    cellInfoAllCells(k, 6, :) = k; % storing cell index
end

%% gathering time-varying intensity for each ROI
parfor j=1:framestop
    im = imread(filename, j);
    % gather pixel intensities and area for each cell
    stats = regionprops(CC, im, 'PixelValues', 'Area');
    % store mean intensity in cellInfoAllCells
    intensities = zeros(size(cellInfoAllCells(:,4,j)));
    for k = 1:CC.NumObjects
        current = stats(k);
        intensities(k) = sum(current.PixelValues) / current.Area;
    end
    
    cellInfoAllCells(:, 4, j) = intensities; % storing raw intensity
end

fprintf('Data Extracted \n')

%% Process data: calculate F-F0/F0
for j = 1:nActiveCells
    % subtract linear trend from mean intensity data
    cellInfoAllCells(j,3,:) = detrend(squeeze(cellInfoAllCells(j,4,:))) + cellInfoAllCells(j,4,1);
    
    % create moving baseline vector
    F = squeeze(cellInfoAllCells(j,4,:));
    window = round(length(F)/100); % moving window is a hundredth of total time
    F = smooth(F, window); % smooth time-varying intensity trace
    F0 = zeros(size(F));
    for k = 1:length(F)
        if k <= window
            F0(k) = prctile(F(1:k+window), 8);
        elseif k >= length(F) - window
            F0(k) = prctile(F(k-window:length(F)), 8);
        else
            F0(k) = prctile(F(k-window:k+window), 8);
        end
    end
    
    % subtract baseline from intensity, then divide by baseline
    cellInfoAllCells(j,3,:) = (squeeze(cellInfoAllCells(j,4,:))-F0)./F0;
    
end

fprintf('Data Processed \n')
cellInfo = cellInfoAllCells; % cellinfo will contain only info for active cells

%% Filter out inactive cells to avoid false positives and measure activity
delete = pi*ones(size(cellInfo(1,3,:))); %marker vector
n = nActiveCells;
for j = 1:n
    signal = smooth(squeeze(cellInfo(j, 3, :)));
    % A spike is a peak whose deltaF/F value is at least 0.25 of the max value in the signal, or 0.05, whichever is higher
    [spikes, ~] = findpeaks(signal, 'MinPeakHeight', max(0.25*max(signal), 0.05));
    
    cellInfo(j, 5, :) = numel(spikes); % storing number of spikes
    cellInfoAllCells(j, 5, :) = numel(spikes); % storing number of spikes
    
    if length(spikes) < 3
        cellInfo(j,3,:) = delete;
        nActiveCells = nActiveCells - 1;
        continue;
    end
    
end

% Delete data from inactive cells in cellinfo (but preserved in cellinfo_allCells)
j = 1;
while size(cellInfo,1) > nActiveCells
    if cellInfo(j,3,:) == delete
        cellInfo(j,:,:) = [];
    else
        j = j + 1;
    end
end
fprintf('Inactive Cells Removed \n')

maxLag = 0; % setting maximum value of lag allowed for calculating correlations
addpath(genpath('/Users/arunmahadevan/Documents/MATLAB/wgPlot'));
loadFile = strrep(filename, '.tif', '_parameters.mat');
load(loadFile, 'cellinfo', 'maxImage', 'mask');

%% Generate adjacency matrix based on cross-correlation of pixel intensities, only for active cells
nActiveCells = size(cellinfo, 1);
adjacencyMatrixAllCorrelations = zeros(nActiveCells, nActiveCells);
cutoff = generateCutoff(cellinfo, 99, maxLag); % using 99th percentile of scrambled data to generate cutoff
fprintf('Cutoff Generated\n')

for j = 1:nActiveCells
    correlations = zeros(length(adjacencyMatrixAllCorrelations), 1);
    parfor k = j:nActiveCells
        current = max(xcov(squeeze(cellinfo(j,3,:)), squeeze(cellinfo(k,3,:)), maxLag, 'coeff'));
        if j == k
            current = 0;
        end
        correlations(k) = current;
    end
    adjacencyMatrixAllCorrelations(j,:) = correlations;
end

functionalAdjacencyMatrix = adjacencyMatrixAllCorrelations;
functionalAdjacencyMatrix(adjacencyMatrixAllCorrelations < cutoff) = 0; % remove indices lower than cutoff
fprintf('Adjacency Matrix Created \n')

%% compute distance between all cell pairs
CellLocations = [squeeze(cellinfo(:, 1, 1)) squeeze(cellinfo(:, 2, 1))]; % using first frame to compute cell distances
CellDistances = squareform(pdist(CellLocations, 'Euclidean'));

saveFile = outputFileName(filename, '_parameters.mat');
save(saveFile, 'functionalAdjacencyMatrix', 'adjacencyMatrixAllCorrelations', 'CellDistances', 'cutoff', '-append');

%% overlay mask on maxImage
f = figure; 
%imshow(segmentOverlay(bw, maxImage));
imshow(maxImage, 'Border', 'Tight'); hold on;
MaskBoundaries = bwboundaries(mask);
for k = 1:length(MaskBoundaries)
    boundary_i = MaskBoundaries{k};
    plot(boundary_i(:, 2), boundary_i(:, 1), 'r');
end

for k = 1:nActiveCells
    text(cellinfo(k,1,1), cellinfo(k,2,1), num2str(cellinfo(k, 6, 1)), 'color', 'w', 'FontSize', 10);
end
saveas(f, strcat(outputFileName(filename, sprintf('-Masks')), '.png'));
close(f);

%% plot graph on maxImage
f = figure; 
%imshow(segmentOverlay(bw, maxImage));
imshow(maxImage); hold on;
MaskBoundaries = bwboundaries(mask);
for k = 1:length(MaskBoundaries)
    boundary_i = MaskBoundaries{k};
    plot(boundary_i(:, 2), boundary_i(:, 1), 'r');
end

for k = 1:nActiveCells
    text(cellinfo(k,1,1), cellinfo(k,2,1), num2str(cellinfo(k, 6, 1)), 'color', 'w', 'FontSize', 6);
end

wgPlot(functionalAdjacencyMatrix, [cellinfo(:,1,1) cellinfo(:,2,1)], 'edgeColorMap', parula); 
saveas(f, strcat(outputFileName(filename, sprintf('-Graph')), '.png'));
close(f);

%% plot graph on immunostain image
if exist(ImmunostainingFileNames{1}, 'file'); blue = imread(ImmunostainingFileNames{1}); else; blue = zeros(size(maxImage)); end
if exist(ImmunostainingFileNames{3}, 'file'); red = imread(ImmunostainingFileNames{3}); else; red = zeros(size(maxImage)); end
if exist(ImmunostainingFileNames{2}, 'file'); green = imread(ImmunostainingFileNames{2}); else; green = zeros(size(maxImage)); end
ImmunostainImage = cat(3, 2*red, zeros(size(green)), blue);

info = imfinfo(filename);
ImmunostainImage = imresize(ImmunostainImage, [info(1).Height, info(1).Width]);

f = figure; set(gcf, 'Color', 'w');
imshow(ImmunostainImage, 'Border', 'Tight');
hold on;
MaskBoundaries = bwboundaries(mask);
hold on;
for k = 1:length(MaskBoundaries)
    boundary_i = MaskBoundaries{k};
    plot(boundary_i(:, 2), boundary_i(:, 1), 'w');
end
%set(findall(gcf,'type','line'), 'LineWidth', 2);

for k = 1:nActiveCells
    text(cellinfo(k,1,1), cellinfo(k,2,1), num2str(cellinfo(k, 6, 1)), 'color', 'w', 'FontSize', 6);
end

wgPlot(functionalAdjacencyMatrix, [cellinfo(:,1,1) cellinfo(:,2,1)], 'edgeColorMap', parula);
saveas(f, strcat(outputFileName(filename, sprintf('-ImmunostainGraph')), '.png'));
close(f);

%% plotting distance v correlation
AvgeDistanceConnectedCells = mean(CellDistances(functionalAdjacencyMatrix~=0));
f = figure;
set(gcf, 'Color', 'w');
if ~isempty(CellDistances)
    plot(CellDistances, adjacencyMatrixAllCorrelations, 'k.', 'MarkerSize', 16);
end

if isnan(cutoff)
    cutoff = 0;
end

%set(gca, 'YLim', [cutoff 1]);
title(sprintf('Average Distance between connected cells = %4.3f pixels', AvgeDistanceConnectedCells));
r = refline(0,cutoff);
r.Color = 'r';
legend(sprintf('Cutoff = %4.3f', cutoff));
xlabel('Distance (pixels)');
ylabel('Correlation Coefficient');
%set(gca, 'FontSize', 30);
saveas(f, strcat(outputFileName(filename, sprintf('-DistanceCorrelation')), '.png'));
close(f);
end