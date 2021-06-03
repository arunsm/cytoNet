
function spatialGraph(outputPath, mask, threshold, adjacencyType)
fprintf('[spatialGraph]Start\n');
% Compute graph based on spatial data
% Adjacency Types
% 1: Centroid distance
% 2: Shared border
globalMetrics = {'', 0};
localMetrics = {'', 0};
adjacencyMatrix = {'', 0};
[globalMetrics{1, 2}, localMetrics{1, 2}, adjacencyMatrix{1, 2}, ...
 maskBoundaries,  cellLocations, globalMetricNames, localMetricNames, ...
 processParams] = ...
  CalculateNetworkMetrics(mask, threshold, adjacencyType);

if ~processParams.processCompleted
    errorReport = makeErrorStruct(['processing for file exceeds computing capacity; spatial graph not created'], 1);
    return;
end

WriteNetworkMetrics(outputPath, adjacencyMatrix, globalMetrics, globalMetricNames, localMetrics, localMetricNames, '_spa');
fprintf('[spatialGraph]Finish\n');
end
