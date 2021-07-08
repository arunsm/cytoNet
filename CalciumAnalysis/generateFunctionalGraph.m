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
mask = imresize(mask, [info(1).Height, info(1).Width]); % ensure mask and images are same dimensions
CC = bwconncomp(mask);
nCells = CC.NumObjects;

if nCells == 0
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
    currentFrame = imread(filePath, j);
    stats = regionprops(CC, currentFrame, 'PixelValues', 'Area');
    intensitiesCurrentFrame = zeros(1, nCells);
    parfor k = 1:nCells
        current = stats(k);
        intensitiesCurrentFrame(k) = sum(current.PixelValues)/current.Area;
    end
    cellInfoAllCells.Fraw(j, :) = intensitiesCurrentFrame;
end

fprintf('Data extracted \n')

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

fprintf('Data processed \n')

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
nActiveCells = sum(cellInfoAllCells.activeCells);
A = zeros(nActiveCells, nActiveCells);
timeSeries = cellInfoAllCells.Fprocessed;

cutoff = generateCutoff(timeSeries, cutoffPercentile, maxLag); % calculate cutoff using scrambled data
cellInfoAllCells.cutoffCorrelation = cutoff;
fprintf('Cutoff generated \n')

for j = 1:nActiveCells
    correlations = zeros(length(A), 1);
    timeSeries1 = timeSeries(:, j);
    parfor k = j:nActiveCells
        timeSeries2 = timeSeries(:, k);
        if j ~= k
            correlations(k) = max(xcov(timeSeries1, timeSeries2, maxLag, 'coeff'));
        end    
    end
    A(j, :) = correlations;
end
functionalAdjacencyMatrixWeighted = A' + A;
functionalAdjacencyMatrixWeighted(1:nActiveCells+1:end) = diag(A);

% create binary adjacency matrix
functionalAdjacencyMatrixBinary = ones(size(functionalAdjacencyMatrixWeighted));
functionalAdjacencyMatrixBinary(functionalAdjacencyMatrixWeighted < cutoff) = 0;
cellInfoAllCells.functionalAdjacencyMatrixBinary = functionalAdjacencyMatrixBinary;
cellInfoAllCells.functionalAdjacencyMatrixWeighted = functionalAdjacencyMatrixWeighted;
fprintf('Adjacency Matrix created \n')
end