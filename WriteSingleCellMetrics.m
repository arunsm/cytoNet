function [] = WriteSingleCellMetrics(OutputPath, MetricNames, Metrics)

nMetrics = numel(MetricNames);
nImages = size(Metrics, 1);

for k = 1:nImages
    
    currentMaskName = Metrics{k, 1};
    
    % opening file for writing local metrics
    WritePath = strcat(OutputPath, filesep, 'SingleCellMetrics_', currentMaskName, '.csv');
    fid = fopen(WritePath, 'w');
    
    % printing metric names
    for i = 1:nMetrics
        fprintf(fid, '%s,', MetricNames{i});
    end
    
    fprintf(fid, '\n');
    formatSpec = '%f,';
    
    % printing image names and metrics
    CurrentMetrics = Metrics{k, 2};
    for i = 1:size(CurrentMetrics, 1)
        fprintf(fid, formatSpec, CurrentMetrics(i, :));
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end
end