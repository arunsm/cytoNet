
% Major and minor axis lengths are scaled by cell area

function [header values] = ComputeNucleusAxisLengthAndRatio(data)

header = {'Nucleus Major Axis Length', 'Nucleus Minor Axis Length', 'Nucleus Major-Minor Axis Length Ratio'};
nucleusmasks = data.nucleusmaskdata;
numnucleusmasks = numel(nucleusmasks);
values = zeros(numnucleusmasks, 3);
for idx_cell = 1:numnucleusmasks
    prop = regionprops(nucleusmasks{idx_cell}, 'Area', 'MajorAxisLength', 'MinorAxisLength');
    largestArea = 0;
    majorAxisLengthLargestArea = 0;
    minorAxisLengthLargestArea = 0;
    for j = 1:numel(prop)
        area = prop(j).Area;
        if area > largestArea
            largestArea = area;
            majorAxisLengthLargestArea = prop(j).MajorAxisLength;
            minorAxisLengthLargestArea = prop(j).MinorAxisLength;
        end
    end
    cellmask = data.cellmaskdata{idx_cell};
    cellprop = regionprops(cellmask, 'Area');
    rootArea = sqrt(cellprop.Area);
    values(idx_cell, 1) = majorAxisLengthLargestArea / rootArea;
    values(idx_cell, 2) = minorAxisLengthLargestArea / rootArea;
    values(idx_cell, 3) = majorAxisLengthLargestArea / minorAxisLengthLargestArea;

end
end
