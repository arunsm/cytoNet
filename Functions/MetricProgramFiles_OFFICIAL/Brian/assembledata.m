function [header data] = assembledata(metricresults, metricfunctions)
% Determine the number of individual metric columns.  Some of the
% metric functions may compute more than one metric (column)
[numrows numcols] = size(metricresults);
nummetrics = 0;
okmetricindex = 0;
for j = 1:numcols
    metric = metricfunctions{j}
    % if there was a problem with the metric function, then skip it
    if metric.skip continue; end
    okmetricindex = j;
    metres = metricresults{1, j}
    if isempty(metres)
        continue;
    else
        nummetrics = nummetrics + size(metres.values, 2);
    end
end

numcells = 0;
if okmetricindex ~= 0
    for i = 1:numrows
        %         display(i)
        %         metres.values
        metres = metricresults{i, okmetricindex};
        if isempty(metres)
            continue;
        else
            numcells = numcells + size(metres.values, 1);
        end
    end
end

% Construct header
header = cell(1, nummetrics);
index = 1;
for j = 1:numcols
    metric = metricfunctions{j};
    if metric.skip continue; end
    metres = metricresults{1, j};
    if isempty(metres)
        continue; 
    end
    localheader = metres.header;
    if strcmp('char', class(localheader))
        header{index} = localheader;
        index = index + 1;
    else
        for k = 1:numel(localheader)
            header{index} = localheader{k};
            index = index + 1;
        end
    end
end

% Construct data table
data = zeros(numcells, nummetrics);
% Copy smaller tables into data table.
% Datarow and datacol keep track of the last filled
% position in data table
datacol = 0;
for j = 1:numcols
    datarow = 0;
    metric = metricfunctions{j};
    if metric.skip continue; end
    for i = 1:numrows
        metres = metricresults{i, j};
        if isempty(metres)
            continue;
        end
        localdata = metres.values;
        [localrows localcols] = size(localdata);
        for ldr = 1:localrows
            for ldc = 1:localcols
                data(datarow + ldr, datacol + ldc) = localdata(ldr, ldc);
            end
        end
        datarow = datarow + localrows;
        
    end
    datacol = datacol + localcols;
end

end




