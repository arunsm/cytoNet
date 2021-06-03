function [header values] = Compute5GLCMMetrics(data)

%fprintf('Compute5GLCMMetrics: Use reduced image for each cell, and not the entire image!!!\n');

metrics = {'Contrast', 'Correlation', 'Energy', 'Entropy from GLCM', 'Homogeneity'};

nummetrics = numel(metrics);
staindata = data.staindata;
numstains = numel(staindata);
stainNames = fieldnames(data.stainkey);
% Limit the number of stain names to the number of stain images
stainNames = stainNames(1:numstains);
header = {};
for i = 1:nummetrics
    h = strcat(metrics{i}, {' '}, stainNames);
    for j = 1:numel(h)
        header{end+1} = h{j};
    end
end
% Add metric labels for rotated cells
rotatedheader = strcat(header, {' '}, 'Rotated');
for i = 1:numel(rotatedheader)
    header{end+1} = rotatedheader{i};
end
% Compute the column distance between a metric for a cell
% and the metric for the rotated cell
offset = numel(header) / 2;


stainkey = data.stainkey;
nucleusIndex = stainkey.nucleus;

maskdata = data.cellmaskdata;
numcells = numel(maskdata);
values = zeros(numcells, numel(header));
for i = 1:numcells
%     fprintf('Cell: %d\n', i);
    mask = maskdata{i};
    
    cc = bwconncomp(mask);
    assert(cc.NumObjects == 1, sprintf('Compute5GLCMMetrics: Image %d in file %s contains %d masks', i, data.cellmaskfilename, cc.NumObjects));
    % Get cell orientation and bounding box of mask and rotated mask
    maskProps = regionprops(mask, 'Orientation', 'BoundingBox');
    bb = maskProps.BoundingBox;
    colstart = ceil(bb(1));
    rowstart = ceil(bb(2));
    numcols = bb(3);
    numrows = bb(4);
    colend = colstart + numcols - 1;
    rowend = rowstart + numrows ;
    orientationAngle = maskProps.Orientation;
    rotatedMask = imrotate(mask, -orientationAngle);
    rotatedMaskProps = regionprops(rotatedMask, 'BoundingBox');
    rotatedbb = rotatedMaskProps.BoundingBox;
    rotatedcolstart = ceil(rotatedbb(1));
    rotatedrowstart = ceil(rotatedbb(2));
    rotatednumcols = rotatedbb(3);
    rotatednumrows = rotatedbb(4);
    rotatedcolend = rotatedcolstart + rotatednumcols - 1;
    rotatedrowend = rotatedrowstart + rotatednumrows ;

    for j = 1:numstains
        stain = staindata{j};
        singlecellstain = mask .* staindata{j};
        
        % For nuclei, compute orientation angle and bounding box from
        % nucleus mask
        if j == nucleusIndex
            nucleusMask = data.nucleusmaskdata{i};
            cc = bwconncomp(nucleusMask);
            % For multinucleated cells, use the largest nucleus
            if cc.NumObjects > 1
                index = largestMaskIndex(cc);
                nucleusMask = logical(zeros(size(nucleusMask)));
                nucleusMask(cc.PixelIdxList{index}) = logical(1);
            end
            nucleusMaskProps = regionprops(nucleusMask, 'Orientation', 'BoundingBox');
            nucleusOrientationAngle = nucleusMaskProps.Orientation;
            rotatedNucleusMask = imrotate(nucleusMask, -nucleusOrientationAngle);
            cc = bwconncomp(rotatedNucleusMask);
            % It is possible for the image rotation to result in more than
            % one connected component.  If so, use the largest.
            if cc.NumObjects > 1
                index = largestMaskIndex(cc);
                rotatedNucleusMask = logical(zeros(size(rotatedNucleusMask)));
                rotatedNucleusMask(cc.PixelIdxList{index}) = 1;
            end
            rotatedNucleusMaskProps = regionprops(rotatedNucleusMask, 'BoundingBox');
            nucleusbb = nucleusMaskProps.BoundingBox;
            nucleuscolstart = ceil(nucleusbb(1));
            nucleusrowstart = ceil(nucleusbb(2));
            nucleusnumcols = nucleusbb(3);
            nucleusnumrows = nucleusbb(4);
            nucleuscolend = nucleuscolstart + nucleusnumcols - 1;
            nucleusrowend = nucleusrowstart + nucleusnumrows - 1;
            rotatednucleusbb = rotatedNucleusMaskProps.BoundingBox;
            rotatednucleuscolstart = ceil(rotatednucleusbb(1));
            rotatednucleusrowstart = ceil(rotatednucleusbb(2));
            rotatednucleusnumcols = rotatednucleusbb(3);
            rotatednucleusnumrows = rotatednucleusbb(4);
            rotatednucleuscolend = rotatednucleuscolstart + rotatednucleusnumcols - 1;
            rotatednucleusrowend = rotatednucleusrowstart + rotatednucleusnumrows - 1;
            % Rotate and isolate stain
            isolatedcellstain = singlecellstain(nucleusrowstart:nucleusrowend, nucleuscolstart:nucleuscolend);
            rotatedsinglecellstain = imrotate(singlecellstain, -nucleusOrientationAngle);
            rotatedisolatedcellstain = rotatedsinglecellstain(rotatednucleusrowstart:rotatednucleusrowend, rotatednucleuscolstart:rotatednucleuscolend);
        else
            % Rotate and isolate stain
            isolatedcellstain = singlecellstain(rowstart:rowend, colstart:colend);
            rotatedsinglecellstain = imrotate(singlecellstain, -orientationAngle);
            rotatedisolatedcellstain = rotatedsinglecellstain(rotatedrowstart:rotatedrowend, rotatedcolstart:rotatedcolend);

        end
        t = statxture_glcm(isolatedcellstain);
        rt = statxture_glcm(rotatedisolatedcellstain);
%        t2 = statxture_glcm(singlecellstain);

        col = j;
        for k = 1:numel(t)
            values(i, col) = t(k);
            values(i, col+offset) = rt(k);
            col = col + numstains;
%            fprintf('%d:  %d  %d\n', k, t(k), t2(k));
        end
    end
end
end


function index = largestMaskIndex(cc)
index = 1;
largestSize = numel(cc.PixelIdxList{1});
for i = 2:cc.NumObjects
    sz = numel(cc.PixelIdxList{i});
    if sz > largestSize
        index = i;
        largestSize = sz;
    end
end
end

%%%
%%% statxture_glcm
%%%
function s = statxture_glcm(f)
offsets = [0 1];
g = graycomatrix(f,'NumLevels',256,'Symmetric', true,'Offset',offsets);
s = struct2cell(ext_graycoprops2(g,'all'));
s = [s{:}];
end
