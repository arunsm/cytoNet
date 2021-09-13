% Function to plot a metric on cell mask image as a heatmap
% INPUT: AdjacencyMatrix (nxn matrix), CellLocations (nx2 matrix), metric
% (nx1 array), Masks (binary image file)
% OUTPUT: figure handle

function f = plotMetricOnImage(cellLocations, masks, metric)

CC = bwconncomp(masks);
nNodes = size(cellLocations, 1);
colors = hot(numel(unique(metric)));
cmap = zeros(nNodes, 3);
L = labelmatrix(CC);
UniqueMetricValues = unique(metric);

for i = 1:numel(UniqueMetricValues)
    index = ismember(metric, UniqueMetricValues(i));
    cmap(index, :) = repmat(colors(i, :), size(cmap(index, :), 1), 1);
end

metricHeatmap = label2rgb(L, cmap);

%f = figure('Visible', 'Off'); 
f = figure;
set(gcf, 'color', 'w');
imshow(metricHeatmap);
hold on;
colormap(hot);
c = colorbar ('northoutside');
c.Ticks = [0 1];
c.TickLabels = [min(metric), max(metric)];
c.FontSize = 18;

end