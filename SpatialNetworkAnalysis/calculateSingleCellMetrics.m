function [SingleCellMetricNames, SingleCellMetrics] = calculateSingleCellMetrics(I, Masks)
    
    I = mat2gray(I);
    SingleCellMetricNames = {'Object Index', 'Centroid(x)', 'Centroid(y)', ...
        'Object Size', 'Circularity', 'Elongation', 'Intensity'};
    
    CC = bwconncomp(Masks);
    nObjects = CC.NumObjects;
    
    ObjectIndex = (1:nObjects)';
    
    c = regionprops(CC, 'Centroid');
    CellLocations = cat(1, c.Centroid);
    
    cellProp = regionprops(CC, 'Area', 'Perimeter');
    ObjectSize = cat(1, cellProp.Area); % Object size is total pixels in object
    ObjectPerimeter = cat(1, cellProp.Perimeter);
    Circularity = 4*pi*ObjectSize./ObjectPerimeter.^2; % circularity = 4*pi*Area/Perimeter^2
    Elongation = ObjectPerimeter./ObjectSize; % elongation = Perimeter/Area
    
    Intensity = zeros(nObjects, 1); % Grayscale intensity of object; number between 0 and 1
    for i = 1:nObjects
        Intensity(i) = mean(I(CC.PixelIdxList{i}));
    end
    
    SingleCellMetrics = [ObjectIndex, CellLocations, ObjectSize, Circularity, Elongation, Intensity];
end