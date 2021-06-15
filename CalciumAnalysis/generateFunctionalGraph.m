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

function [errorReport, cellInfo, cellInfoAllCells, mask, maxImage, functionalAdjacencyMatrix] = generateFunctionalGraph(filename, mask)

%% Setting up
nSpikesThreshold = 2;
info = imfinfo(filename);
maxImage = findMaxImage(filename); % extract maximum intensity projection image

% if mask not provided as input, run default segmentation
%if nargin < 2
    mask = segmentWatershed(maxImage);
%end

totalFrames = size(info, 1);
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

%% gathering information on ROIs
nActiveCells = nCells; % initializing number of active cells to all cells
cellInfoAllCells = struct;
cellInfoAllCells.nCells = nCells;
stats = regionprops(CC, 'Centroid');

% centroid locations of all cells
cellInfoAllCells.cellLocations = zeros(nCells, 2);
for k = 1:nCells
    current = stats(k);
    location = current.Centroid;
    cellInfoAllCells.cellLocations(k, 1) = location(1); % storing x location
    cellInfoAllCells.cellLocations(k, 2) = location(2); % storing y location
end

% store unprocessed time-varying intensity for each ROI
cellInfoAllCells.Fraw = zeros(totalFrames, nCells);
for j = 1:totalFrames
    currentFrame = imread(filename, j);
    stats = regionprops(CC, currentFrame, 'PixelValues', 'Area');
    intensitiesCurrentFrame = zeros(1, nCells);
    parfor k = 1:nCells
        current = stats(k);
        intensitiesCurrentFrame(k) = sum(current.PixelValues)/current.Area;
    end
    cellInfoAllCells.Fraw(j, :) = intensitiesCurrentFrame;
end

fprintf('Data Extracted \n')

%% Process data: calculate F-F0/F0

cellInfoAllCells.Fprocessed = zeros(totalFrames, nCells);
cellInfoAllCells.F0 = zeros(totalFrames, nCells);
Fraw = cellInfoAllCells.Fraw;

for j = 1:nCells % loop over all cells
    currentF = Fraw(:, j);
    
    % subtract linear trend from mean intensity data
    currentF = detrend(currentF) + currentF(1);
    
    % create moving baseline vector
    window = round(totalFrames/100); % moving window is a hundredth of total time
    currentF = smooth(currentF, window); % smooth time-varying intensity trace
    F0 = zeros(size(currentF));
    parfor k = 1:totalFrames
        if k <= window
            F0(k) = prctile(currentF(1:k+window), 8);
        elseif k >= length(currentF) - window
            F0(k) = prctile(currentF(k-window:length(currentF)), 8);
        else
            F0(k) = prctile(currentF(k-window:k+window), 8);
        end
    end
    
    cellInfoAllCells.F0(:, j) = F0; % store baseline intensity
    cellInfoAllCells.Fprocessed(:, j) = (currentF - F0)./F0; % calculate delta F / F and store as Fprocessed
end

fprintf('Data Processed \n')

%% Filter out inactive cells to avoid false positives and measure activity
cellInfoAllCells.activeCells = zeros(1, nCells); % marker vector for active cells
cellInfoAllCells.nSpikes = zeros(1, nCells); % storing number of spikes for each cell
for j = 1:nCells
    currentSignal = cellInfoAllCells.Fprocessed(:, j);
    
    % A spike is a peak whose deltaF/F value is at least 0.25 of the max value in the signal, or 0.05, whichever is higher
    [spikes, ~] = findpeaks(currentSignal, 'MinPeakHeight', max(0.25*max(currentSignal), 0.05));
    
    cellInfoAllCells.nSpikes(j) = numel(spikes); % storing number of spikes
    
    % setting marker vector to 1 if number of spikes exceeds threshold
    if numel(spikes) > nSpikesThreshold
        cellInfoAllCells.activeCells(j) = 1;
    end
end

%% Generate adjacency matrix based on cross-correlation of pixel intensities, only for active cells
maxLag = 0; % setting maximum value of lag allowed for calculating correlations
addpath(genpath('/Users/arunmahadevan/Documents/MATLAB/wgPlot'));
loadFile = strrep(filename, '.tif', '_parameters.mat');
load(loadFile, 'cellinfo', 'maxImage', 'mask');

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