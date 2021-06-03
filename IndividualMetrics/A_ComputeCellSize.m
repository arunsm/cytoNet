function [header values] = ComputeCellSize(data)

fprintf('ComputeCellSize:  Needs to be scaled by magnification level!!!\n');

header = 'Cell Size';
masks = data.cellmaskdata;
nummasks = numel(masks);
values = zeros(nummasks, 1);
for i = 1:nummasks
    cellProp = regionprops(masks{i}, 'Area');
    values(i) = cellProp.Area;
    end
end