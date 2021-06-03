function [header values] = ComputeNucleusMajorAxisLength(data)
% fprintf('ComputeNucleusMajorAxisLength: scale this by root of cell area!!!\n');

header = 'Nucleus Major Axis Length';
nucleusmasks = data.nucleusmaskdata;
numnucleusmasks = numel(nucleusmasks);
values = zeros(numnucleusmasks, 1);
for idx_cell = 1:numnucleusmasks
    prop = regionprops(nucleusmasks{idx_cell}, 'Area', 'MajorAxisLength');
    largestArea = 0;
    majorAxisLengthLargestArea = 0;
    for j = 1:numel(prop)
        area = prop(j).Area;
        if area > largestArea
            largestArea = area;
            majorAxisLengthLargestArea = prop(j).MajorAxisLength;
        end
    end
    cellmask = data.cellmaskdata{idx_cell};
    cellprop = regionprops(cellmask, 'Area');
    values(idx_cell) = majorAxisLengthLargestArea / sqrt(cellprop.Area);
end
end