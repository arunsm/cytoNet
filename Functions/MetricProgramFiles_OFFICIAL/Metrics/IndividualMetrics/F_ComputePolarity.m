function [header values] = ComputePolarity(data)
%fprintf('ComputePolarity: Scale polarity distance by square root of cell area\n');

cellMask = data.cellmaskdata;
numOfCell = numel(cellMask);
stainImage = data.staindata;
numOfStain = numel(stainImage);

stainNames = fieldnames(data.stainkey);
% Putting the curly braces {} around the polarity string prevents its
% trailing space from being removed by strcat
header = strcat({'Polarity '}, stainNames(1:numOfStain))';
values = zeros(numOfCell, numOfStain);
for idx_cell = 1:numOfCell
    cellProp = regionprops(cellMask{idx_cell},'PixelIdxList', 'PixelList','Centroid', 'Area');
    rootarea = sqrt(cellProp.Area);
    idxPixel = cellProp.PixelIdxList;
    xCell  = cellProp.Centroid(1);
    yCell  = cellProp.Centroid(2);
    xStain = cellProp.PixelList(:,1);
    yStain = cellProp.PixelList(:,2);
    for idx_stain = 1:numOfStain
        stainImageOneCell = cellMask{idx_cell}.*stainImage{idx_stain};
        pixelValues = stainImageOneCell(idxPixel);
        sumPixelValues = sum(pixelValues);  
        xBarStain = sum(xStain.* pixelValues) / sumPixelValues;
        yBarStain = sum(yStain.* pixelValues) / sumPixelValues;
        values(idx_cell,idx_stain) = sqrt((xBarStain-xCell).^2+(yBarStain-yCell).^2) / rootarea;
    end
end
end