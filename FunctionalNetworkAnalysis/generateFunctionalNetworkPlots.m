function errorReport = generateFunctionalNetworkPlots(cellInfoAllCells, fileName, outputDirName)

errorReport = [];
fileNameBase = createBaseName(fileName);

%% correlation between distance and functional edge strength
cellLocationsConnectedCells = cellInfoAllCells.cellLocations(find(cellInfoAllCells.activeCells), :);
cellDistancesConnectedCells = squareform(pdist(cellLocationsConnectedCells, 'Euclidean'));
avgeDistanceConnectedCells = mean(cellDistancesConnectedCells(cellInfoAllCells.adjacencyMatrixWeighted~=0));

f = figure('Visible', 'Off');
set(gcf, 'Color', 'w');
if ~isempty(cellDistancesConnectedCells)
    plot(cellDistancesConnectedCells, cellInfoAllCells.adjacencyMatrixWeighted, 'k.', 'MarkerSize', 16);
end

if isnan(cellInfoAllCells.cutoffCorrelation)
    cellInfoAllCells.cutoffCorrelation = 0;
end

title(sprintf('Average Distance between connected cells = %4.3f pixels', avgeDistanceConnectedCells));
if ~isnan(cellInfoAllCells.cutoffCorrelation)
    r = refline(0, cellInfoAllCells.cutoffCorrelation);
    r.Color = 'r';
end

xlabel('Distance (pixels)');
ylabel('Functional Edge Strength');
set(gca, 'FontSize', 20);
savePath = strcat(outputDirName, filesep, fileNameBase, '-DistanceCorrelation.svg');
saveas(f, savePath);
close(f);

%% plot weighted graph on maxImage (only correlations above cutoff)
f = figure('Visible', 'Off');
set(gcf, 'Color', 'w');
imshow(cellInfoAllCells.maxImage); hold on;
maskBoundaries = bwboundaries(cellInfoAllCells.mask);
for k = 1:length(maskBoundaries)
    currentBoundary = maskBoundaries{k};
    plot(currentBoundary(:, 2), currentBoundary(:, 1), 'r');
end

for k = 1:cellInfoAllCells.nNodes
    text(cellInfoAllCells.cellLocations(k, 1), cellInfoAllCells.cellLocations(k, 2), num2str(k), 'color', 'r', 'FontSize', 8);
end

A = cellInfoAllCells.adjacencyMatrixWeighted;
A(cellInfoAllCells.adjacencyMatrixWeighted < cellInfoAllCells.cutoffCorrelation) = 0;
wgPlot(A, cellInfoAllCells.cellLocations(find(cellInfoAllCells.activeCells), :), ...
    'edgeColorMap', parula);
set(findall(gcf,'type','line'), 'LineWidth', 1);
savePath = strcat(outputDirName, filesep, fileNameBase, '-FunctionalGraph.svg');
saveas(f, savePath);
close(f);

%% heatmap of number of spikes overlaid on mask

f = plotMetricOnImage(cellInfoAllCells.cellLocations, cellInfoAllCells.mask, cellInfoAllCells.nSpikes);
hold on;

for k = 1:cellInfoAllCells.nNodes
    text(cellInfoAllCells.cellLocations(k, 1), cellInfoAllCells.cellLocations(k, 2), num2str(k), 'color', 'k', 'FontSize', 8);
end

savePath = strcat(outputDirName, filesep, fileNameBase, '-NumberSpikes.svg');
saveas(f, savePath);
close(f);

end