function [header values] = ComputeCircularityElongation(data)
header = {'Circularity', 'Elongation'};
cellmasks = data.cellmaskdata;
nummasks = numel(cellmasks);
values = zeros(nummasks, 2);
for i = 1:nummasks
    cellProp = regionprops(cellmasks{i}, 'Area', 'Perimeter');
    values(i, 1) = 4 * pi * cellProp.Area / (cellProp.Perimeter)^2;
    values(i, 2) = cellProp.Perimeter/cellProp.Area;
end
end