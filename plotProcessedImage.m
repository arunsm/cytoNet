
function [] = plotProcessedImage(outputPath, cellInfoAllCells)

fileName = cellInfoAllCells.fileName;
adjacencyMatrixBinary = cellInfoAllCells.adjacencyMatrixBinary;
mask = cellInfoAllCells.mask;
cellLocations = cellInfoAllCells.cellLocations;
graphTypeTag = cellInfoAllCells.graphTypeTag;

writePath = strcat(outputPath, filesep, fileName, '-PROCESSED-', graphTypeTag, '.tif');

% compute cell locations and ROI boundaries from mask
maskBoundaries = bwboundaries(mask);

f = figure('Visible', 'Off');
imshow(mask, 'Border', 'Tight'); hold on;

if ~isempty(adjacencyMatrixBinary)
    gplot(adjacencyMatrixBinary, cellLocations, 'y');
end

for k = 1:length(maskBoundaries)
    currentBoundary = maskBoundaries{k};
    plot(currentBoundary(:, 2), currentBoundary(:, 1), 'r');
    text(cellLocations(k, 1), cellLocations(k, 2), num2str(k), 'FontSize', 8, 'Color', 'red');
end

saveas(f, writePath);
close(f);

end