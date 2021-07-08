% Function to write network metrics to file
% INPUT: OutputPath - directory path for output
% AdjacencyMatrices, GlobalMetrics, LocalMetrics - {nx2} cell arrays 
% where n is the number of images and the second column contains the
% adjacency matrix, global metrics array, and local metrics array for
% corresponding image
% GlobalMetricNames, LocalMetricNames - cell arrays containing names of
% metrics
% graphTypeTag - string containing type of graph - 'spa' for spatial, 'fun'
% for functional

function [] = WriteNetworkMetrics(OutputPath, AdjacencyMatrices, GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames, graphTypeTag)

if nargin < 7
    graphTypeTag = '_spa';
end

nGlobalMetrics = numel(GlobalMetricNames);
nLocalMetrics = numel(LocalMetricNames);
nImages = size(GlobalMetrics, 1);

%% Print Global Metrics

% open file for writing global metrics
WritePath = strcat(OutputPath, filesep, 'GlobalMetrics', graphTypeTag, '.csv');
fid = fopen(WritePath, 'w');

% print metric names
fprintf(fid, ',');
for i = 1:nGlobalMetrics
    fprintf(fid, '%s,', GlobalMetricNames{i});
end

fprintf(fid, '\n');
formatSpec = '%f,';

% print image names and metrics
for i = 1:nImages
    currentMaskName = GlobalMetrics{i, 1};
    GlobalMetricsCurrent = GlobalMetrics{i, 2};
    imageName = sprintf('%s', currentMaskName);
    fprintf(fid, '%s,', imageName);
    fprintf(fid, formatSpec, GlobalMetricsCurrent(:));
    fprintf(fid, '\n');
end

fclose(fid);

%% Print Local Metrics

for k = 1:nImages
    
    currentMaskName = LocalMetrics{k, 1};
    
    % opening file for writing local metrics
    WritePath = strcat(OutputPath, filesep, 'LocalMetrics_', currentMaskName, graphTypeTag, '.csv');
    fid = fopen(WritePath, 'w');
    
    % printing metric names
    for i = 1:nLocalMetrics
        fprintf(fid, '%s,', LocalMetricNames{i});
    end
    
    fprintf(fid, '\n');
    formatSpec = '%f,';
    
    % printing image names and metrics
    LocalMetricsCurrent = LocalMetrics{k, 2};
    for i = 1:size(LocalMetricsCurrent, 1)
        fprintf(fid, formatSpec, LocalMetricsCurrent(i, :));
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end

%% Print Adjcacency Matrix

for k = 1:nImages
    
    currentMaskName = AdjacencyMatrices{k, 1};
    
    % opening file for writing local metrics
    WritePath = strcat(OutputPath, filesep, 'AdjacencyMatrix_', currentMaskName, graphTypeTag, '.csv');
    csvwrite(WritePath, AdjacencyMatrices{k, 2});
end
