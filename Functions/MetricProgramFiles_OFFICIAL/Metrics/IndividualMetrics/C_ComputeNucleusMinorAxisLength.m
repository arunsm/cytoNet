function [header values] = ComputeNucleusMinorAxisLength(data)
% fprintf('ComputeNucleusMinorAxisLength: scale this by root of cell area!!!\n');

header = 'Nucleus Minor Axis Length';
nucleusmasks = data.nucleusmaskdata;
numnucleusmasks = numel(nucleusmasks);
values = zeros(numnucleusmasks, 1);
for idx_cell = 1:numnucleusmasks
    prop = regionprops(nucleusmasks{idx_cell}, 'Area', 'MinorAxisLength');
    largestArea = 0;
    minorAxisLengthLargestArea = 0;
    for j = 1:numel(prop)
        area = prop(j).Area;
        if area > largestArea
            largestArea = area;
            minorAxisLengthLargestArea = prop(j).MinorAxisLength;
        end
    end
    cellmask = data.cellmaskdata{idx_cell};
    cellprop = regionprops(cellmask, 'Area');
    values(idx_cell) = minorAxisLengthLargestArea / sqrt(cellprop.Area);  
    end
end