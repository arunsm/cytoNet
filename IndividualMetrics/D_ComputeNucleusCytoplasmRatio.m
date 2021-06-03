function [header values] = ComputeNucleusCytoplasmRatio(data)
%fprintf('ComputeNucleusCytoplasmRatio: What should happen to multi-nucleus cells?\n');
%fprintf('ComputeNucleusCytoplasmRatio: Consider using the cell mask instead of the cytoplasm mask\n');
header = 'Nucleus Cytoplasm Ratio';
nuMask = data.nucleusmaskdata;
% cyMask = data.cytomaskdata;
cellMask = data.cellmaskdata;
assert(numel(nuMask) == numel(cellMask), sprintf('ComputeNucleusCytoplasmRatio: files %s and %s have different number of images',...
    data.nucleusmaskfilename, data.cellmaskfilename));
values = zeros(length(nuMask), 1);
for idx_cell = 1:length(nuMask)
    nuProp = regionprops(nuMask{idx_cell},'Area');
    totalNucleusArea = 0;
    for i = 1:numel(nuProp)
        totalNucleusArea = totalNucleusArea + nuProp(i).Area;
    end
%     cyProp = regionprops(cyMask{idx_cell}, 'Area');
%     totalCytoArea = 0;
%     for i = 1:numel(cyProp)
%         totalCytoArea = totalCytoArea + cyProp(i).Area;
%     end
    cellProp = regionprops(cellMask{idx_cell}, 'Area');
    assert(numel(cellProp) == 1, sprintf('ComputeNucleusCytoplasmRatio: Image %d of file %s has more than on cell',...
        idx_cell, data.cellmaskfilename));
    cellArea = cellProp.Area;
    values(idx_cell) = totalNucleusArea / (cellArea);        %(cellArea - totalNucleusArea);
%     values(idx_cell) = totalNucleusArea / totalCytoArea;
end
end
