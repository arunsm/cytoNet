% Function to extract calcium time series data from an image stack and
% compute a functional graph based on correlations among time series
%
% INPUT: filePath - stacked tif image file you wish to generate a graph for; 
% mask (optional) - binary mask defining ROIs for time series extraction; 
% if no mask is provided, the default segmentation program is run on 
% the maximum intensity projection image
%
% OUTPUT: cellInfoAllCells - MATLAB structure containing calcium image
% analysis parameters

function [errorReport, cellInfoAllCells] = generateFunctionalGraph(filePath, mask)

errorReport = [];
warning('off','all')

%% Setting up
nSpikesThreshold = 2; % minimum number of spikes to determine an active cell
cutoffPercentile = 99; % percentile of scrambled correlations to be used for cutoff
maxLag = 0; %  maximum value of lag allowed for calculating correlations

info = imfinfo(filePath);
maxImage = findMaxImage(filePath); % extract maximum intensity projection image

% if mask not provided as input, run default segmentation
if nargin < 2
    mask = segmentWatershed(maxImage);
end

totalFrames = size(info, 1);
CC = bwconncomp(mask);
nNodes = CC.NumObjects;

if nNodes == 0
    fn = getFileName(filePath);
    errorReport = makeErrorStruct(['no cells found in file ', fn], 1);
    cellInfo = [];
    cellInfoAllCells = [];
    maxImage = [];
    A = [];
    return;
end

%% gathering information on ROIs
cellInfoAllCells = struct;
cellInfoAllCells.mask = mask;
cellInfoAllCells.maxImage = maxImage;
cellInfoAllCells.nNodes = nNodes;
cellInfoAllCells.graphTypeTag = 'functional';

stats = regionprops(CC, 'Centroid');

% centroid locations of all cells
cellInfoAllCells.cellLocations = zeros(nNodes, 2);
for k = 1:nNodes
    current = stats(k);
    location = current.Centroid;
    cellInfoAllCells.cellLocations(k, 1) = location(1); % storing x location
    cellInfoAllCells.cellLocations(k, 2) = location(2); % storing y location
end

% store unprocessed time-varying intensity for each ROI
cellInfoAllCells.Fraw = zeros(totalFrames, nNodes);
for j = 1:totalFrames
    currentFrame = imread(filePath, j);
    stats = regionprops(CC, currentFrame, 'PixelValues', 'Area');
    intensitiesCurrentFrame = zeros(1, nNodes);
    parfor k = 1:nNodes
        current = stats(k);
        intensitiesCurrentFrame(k) = sum(current.PixelValues)/current.Area;
    end
    cellInfoAllCells.Fraw(j, :) = intensitiesCurrentFrame;
end

fprintf('Data extracted \n')

%% Process data: calculate F-F0/F0

cellInfoAllCells.Fprocessed = zeros(totalFrames, nNodes);
cellInfoAllCells.F0 = zeros(totalFrames, nNodes);
Fraw = cellInfoAllCells.Fraw;

for j = 1:nNodes % loop over all cells
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
    cellInfoAllCells.Fprocessed(:, j) = smooth((currentF - F0)./F0); % calculate delta F / F, smooth data, and store as Fprocessed
end

fprintf('Data processed \n')

%% Filter out inactive cells to avoid false positives and measure activity
cellInfoAllCells.activeCells = zeros(1, nNodes); % marker vector for active cells
cellInfoAllCells.nSpikes = zeros(1, nNodes); % storing number of spikes for each cell
for j = 1:nNodes
    currentSignal = cellInfoAllCells.Fprocessed(:, j);
    
    % A spike is a peak whose deltaF/F value is at least 0.25 of the max value in the signal, or 0.05, whichever is higher
    [spikes, ~] = findpeaks(currentSignal, 'MinPeakHeight', max(0.25*max(currentSignal), 0.05), 'MinPeakProminence', max(0.25*max(currentSignal), 0.05));
    
    cellInfoAllCells.nSpikes(j) = numel(spikes); % storing number of spikes
    
    % setting marker vector to 1 if number of spikes exceeds threshold
    if numel(spikes) > nSpikesThreshold
        cellInfoAllCells.activeCells(j) = 1;
    end
end

%% Generate adjacency matrix based on cross-correlation of pixel intensities, only for active cells
nActiveCells = sum(cellInfoAllCells.activeCells);
A = zeros(nActiveCells, nActiveCells);
timeSeries = cellInfoAllCells.Fprocessed(:, find(cellInfoAllCells.activeCells));

cutoff = generateCutoff(timeSeries, cutoffPercentile, maxLag); % calculate cutoff using scrambled data
cellInfoAllCells.cutoffCorrelation = cutoff;
fprintf('Cutoff generated \n')

% for j = 1:nActiveCells
%     correlations = zeros(length(A), 1);
%     timeSeries1 = timeSeries(:, j);
%     parfor k = 1:nActiveCells
%         timeSeries2 = timeSeries(:, k);
%         if j ~= k
%             correlations(k) = max(xcov(timeSeries1, timeSeries2, maxLag, 'coeff'));
%         end
%     end
%     A(j, :) = correlations;
% end
% functionalAdjacencyMatrixWeighted = A' + A;
% functionalAdjacencyMatrixWeighted(1:nActiveCells+1:end) = diag(A);

A = corr(timeSeries);
adjacencyMatrixWeighted = A - diag(diag(A)); % set diagonal elements to zero

% create binary adjacency matrix
adjacencyMatrixBinary = ones(size(adjacencyMatrixWeighted));
adjacencyMatrixBinary(adjacencyMatrixWeighted < cutoff) = 0;
cellInfoAllCells.functionalAdjacencyMatrixBinary = adjacencyMatrixBinary;
cellInfoAllCells.functionalAdjacencyMatrixWeighted = adjacencyMatrixWeighted;
fprintf('Adjacency Matrix created \n')
end