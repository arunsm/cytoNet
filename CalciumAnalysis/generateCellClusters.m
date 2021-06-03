function [idx] = generateCellClusters(cellinfo, cellinfo_extended, ... 
    radius_um, minmembers, pixelwidth, image, fileName, outputDirName)
% The brief is : 
% 4) Define ?cluster cells? (cells create a group of activation) according 
% to specific criteria like within 50um 3 or more cells activation;

% ZM : Standard clustering methods work differently -- either ID the number
% of clusters you want (as in k-means) or use neighboring. Figuring out
% clusters in which all members are within a certain radius is possible,
% but computationally expensive. I think you'd have to determine all the
% members based on radius from each cell with that cell at the center 
% (easier), or use each cell as a center for a circle from which you draw
% potential cluster regions. Will have to see if that's what he's looking
% for.

radiusCluster_px = radius_um / pixelwidth; %um converted to px at pixelwidth / um.
celldata = cellinfo_extended.centroid; 
idx = dbscan(celldata, radiusCluster_px, minmembers);
figure('visible', 'off');
imshow(image,'InitialMagnification',200); 
hold on; 
gscatter(celldata(:,1), celldata(:,2), idx); 
F = getframe; 
[filepath,name,ext] = fileparts(fileName);
imwrite(F.cdata, sprintf('%s%s%s_cellclusters.tif', outputDirName, filesep, name)); 
close(figure); 
end

