function errorReport = calciumEngine(filePath, outputDirName)

fileName = getFileName(filePath);    
fileNameBase = createBaseName(fileName);

fprintf('Generating functional graph... \n')
[errorReport, cellInfoAllCells] = generateFunctionalGraph(filePath);

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
    calculateNetworkMetrics(cellInfoAllCells.functionalAdjacencyMatrixBinary);

if ~processParams.processCompleted
    errorReport = makeErrorStruct(['processing for file ', fileName, ' exceeds computing capacity; graph not created'], 1);
end

% write metrics to file
WriteNetworkMetrics(outputDirName, AdjacencyMatrices, ...
    GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames, '_fun');

% create spatial graph
adjacencyType = 1;
threshold = 1.5;

if isempty(errorReport)
    spatialGraph(outputDirName, bw, threshold, adjacencyType);
end

end