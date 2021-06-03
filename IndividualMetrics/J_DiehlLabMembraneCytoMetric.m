function [header table] = membraneCytoMetric(data)
membranePixelThickness = [3, 5, 10];

numCells = numel(data.cellmaskdata);
numStains = numel(data.staindata);

key = data.stainkey;
nucleusIndex = key.nucleus;
stainNames = fieldnames(key);
% Sort stainNames by their associated values in key
stainIndex = zeros(numel(stainNames), 1);
for i = 1:numel(stainNames)
    stainIndex(i) = getfield(key, stainNames{i});
end
[~, I] = sort(stainIndex);
orderedStainNames = cell(numel(stainNames), 1);
for i= 1:numel(stainNames)
    orderedStainNames{i} = stainNames{I(i)};
end

headerByStainByMembraneSize = {'Membrane Intensity Mean', 'Cytoplasm Intensity Mean',...
    'Membrane Intensity Summed', 'Cytoplasm Intensity Summed',...
    'Membrane Intensity Median', 'Cytoplasm Intensity Median',...
    'Membrane Intensity Std Dev', 'Cytoplasm Intensity Std Dev'};

headerByMembraneSize = {'Membrane Area Total Area Ratio'};

numColumns = (numel(headerByMembraneSize) + (numel(headerByStainByMembraneSize) * (numStains - 1))) * numel(membranePixelThickness);
header = cell(1, numColumns);

idx = 1;
for p = membranePixelThickness
    membraneTag = sprintf('MembThckns%d', p);
    membraneSizeColumnNames = strcat(membraneTag, {' '}, headerByMembraneSize);
    last = idx + numel(membraneSizeColumnNames) - 1;
    header(1, idx:last) = membraneSizeColumnNames;
    idx = last + 1;
    for s = 1:numel(orderedStainNames)
        if s == key.nucleus continue; end
        stainColumnNames = strcat(membraneTag, {' '}, orderedStainNames{s}, {' '},  headerByStainByMembraneSize);
        last = idx + numel(stainColumnNames) - 1;
        header(idx:last) = stainColumnNames;
        idx = last + 1;
    end
end

table = zeros(numCells, numColumns);
for c = 1:numel(data.cellmaskdata)
    idx = 1;
    mask = data.cellmaskdata{c};
    for p = membranePixelThickness
        cytoplasmMask = imerode(mask, strel('disk', p, 0));
        membraneMask = mask & ~cytoplasmMask;

        cytoplasmArea = sum(double(cytoplasmMask(:)));
        membraneArea = sum(double(membraneMask(:)));

        membraneRatio = membraneArea / (membraneArea + cytoplasmArea);

        table(c, idx) = membraneRatio;
        idx = idx + 1;

        for s = 1:numel(data.staindata)
            if s == key.nucleus continue; end

%            cytoplasmIntensity = data.staindata{s} .* double(cytoplasmMask);
%            membraneIntensity = data.staindata{s} .* double(membraneMask);

            % Create a column vector containing only those pixels values of
            % cytoplasm stain intensity and membrane stain intensity.
            cytoplasmIntensity = data.staindata{s}(cytoplasmMask);
            membraneIntensity = data.staindata{s}(membraneMask);

% fprintf('number of cytoplasmIntensity pixels: %d\n', numel(cytoplasmIntensity));
% fprintf('number of membraneIntensity pixels: %d\n', numel(membraneIntensity));

            summedCytoplasmIntensity = sum(cytoplasmIntensity(:));
            summedMembraneIntensity = sum(membraneIntensity(:));

            meanCytoplasmIntensity = summedCytoplasmIntensity / cytoplasmArea;
            meanMembraneIntensity = summedMembraneIntensity / membraneArea;

            medianCytoplasmIntensity = median(cytoplasmIntensity);
            medianMembraneIntensity = median(membraneIntensity);
            
            sdevCytoplasmIntensity = std(cytoplasmIntensity(:));
            sdevMembraneIntensity = std(membraneIntensity(:));

            last = idx + numel(headerByStainByMembraneSize) - 1;
            table(c, idx:last) = [meanMembraneIntensity, meanCytoplasmIntensity,...
                summedMembraneIntensity, summedCytoplasmIntensity,...
                medianMembraneIntensity, medianCytoplasmIntensity,...
                sdevMembraneIntensity, sdevCytoplasmIntensity];

            idx = last + 1;
        end
    end
end

end
