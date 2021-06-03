
function errorReport = writeGraphProperties(OutputPath, adj, CurrentMask, fileName)
errorReport = [];
flag = false;

%% starting parallel pool
%myCluster = parcluster('local');
%myCluster.JobStorageLocation = pwd;
%poolObj = parpool(myCluster);

%% Loading parameters and initializing variables

GlobalMetrics = {}; % global metric cell array
LocalMetrics = {}; % local metric cell array
AdjacencyMatrices = {}; % cell array to hold adjacency matrices

[numRows, numCols] = size(CurrentMask);
%    imageDimensions = [info.Width info.Height];
imageDimensions = [numCols, numRows];
        
flag = true; % this means we have usable input
            
fprintf('\t creating graph representation for image\n');

currentMaskName = 'MASK';
imwrite(CurrentMask, [OutputPath, filesep, currentMaskName, '.tif'], 'tif', 'Compression', 'none');
GlobalMetrics{1, 1} = currentMaskName;
LocalMetrics{1, 1} = currentMaskName;
AdjacencyMatrices{1, 1} = currentMaskName;
        
% calling CalculateNetworkMetrics to compute graph-based metrics
[GlobalMetrics{1, 2}, LocalMetrics{1, 2}, AdjacencyMatrices{1, 2}, MaskBoundaries, ...
    CellLocations, GlobalMetricNames, LocalMetricNames, processParams]...
    = graphMetrics(CurrentMask, adj);
        
if processParams.processCompleted
    % call function to plot processed image figure
    prImagePath = strcat(OutputPath, filesep, currentMaskName, '-PROCESSED.tif');
    plotProcessedImage(prImagePath, CurrentMask, AdjacencyMatrices{1, 2}, CellLocations, MaskBoundaries);
else
    errorReport = makeErrorStruct(['processing for file ', fileName, ' exceeds computing capacity; graph not created'], 1);
end
        
%% writing metrics to file

fprintf('\t\t\twriting metrics to file\n');

if flag
    WriteNetworkMetrics(OutputPath, AdjacencyMatrices, GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames, '_fun');
end

%% cleaning up
%delete(poolObj);

end
