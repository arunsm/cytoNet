
function [] = plotProcessedImage(prImagePath, Masks, AdjacencyMatrix, CellLocations, MaskBoundaries)

f = figure('Visible', 'Off');
imshow(Masks, 'Border', 'Tight'); hold on;

if ~isempty(AdjacencyMatrix)
    gplot(AdjacencyMatrix, CellLocations, 'y');
end

for k = 1:length(MaskBoundaries)
    boundary_i = MaskBoundaries{k};
    plot(boundary_i(:, 2), boundary_i(:, 1), 'r');
    text(CellLocations(k, 1), CellLocations(k, 2), num2str(k), 'FontSize', 8, 'Color', 'red');
end

saveas(f, prImagePath);
close(f);

end