% Script calling functions to generate functional graph, generate
% functional network plots, compute graph-based metrics, and write metrics
% to file
% INPUT: filePath - directory path to input file, assumed to be a stacked
% .tif file containing calcium imaging data; outputDirName - directory path
% to folder for output
% DEPENDENCIES:
% Mike Wu (2021). wgPlot - Weighted Graph Plot (a better version of gplot)
% https://www.mathworks.com/matlabcentral/fileexchange/24035-wgplot-weighted-graph-plot-a-better-version-of-gplot
% MATLAB Central File Exchange

function errorReport = calciumEngine(filePath, outputDirName, optionalMask)

info = imfinfo(filePath);
fileName = getFileName(filePath);
fileNameBase = createBaseName(fileName);

fprintf('Generating functional graph... \n')

if nargin < 3
    [errorReport, cellInfoAllCells] = generateFunctionalGraph(filePath);
else
    if size(optionalMask, 1) ~= info(1).Height || size(optionalMask, 2) ~= info(1).Width
        errorReport = makeErrorStruct(['dimensions of mask and image file do not match', fileName], 1);
        return;
    end
    [errorReport, cellInfoAllCells] = generateFunctionalGraph(filePath, optionalMask);
end

if ~isempty(errorReport)
    return;
end

if sum(cellInfoAllCells.activeCells) == 0
    errorReport = makeErrorStruct(['no active cells found in file ', fileName], 1);
    return;
end

generateFunctionalNetworkPlots(cellInfoAllCells, fileName, outputDirName);

GlobalMetrics = {}; LocalMetrics = {};
GlobalMetrics{1, 1} = fileNameBase; LocalMetrics{1, 1} = fileNameBase;
AdjacencyMatrices = {}; AdjacencyMatrices{1, 1} = fileNameBase; AdjacencyMatrices{1, 2} = cellInfoAllCells.functionalAdjacencyMatrixBinary;

% compute graph-based metrics
[GlobalMetrics{1, 2}, LocalMetrics{1, 2}, GlobalMetricNames, LocalMetricNames, processParams] = ...
    calculateNetworkMetrics(cellInfoAllCells.mask, cellInfoAllCells.functionalAdjacencyMatrixBinary);

if ~processParams.processCompleted
    errorReport = makeErrorStruct(['processing for file ', fileName, ' exceeds computing capacity; graph not created'], 1);
end

% write metrics to file
WriteNetworkMetrics(outputDirName, AdjacencyMatrices, ...
    GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames, '_fun');

% write cellInfoAllCells to file
savePath = strcat(outputDirName, filesep, fileNameBase, '-cellInfoAllCells.m');
save(savePath, 'cellInfoAllCells');

% create spatial graph?
end