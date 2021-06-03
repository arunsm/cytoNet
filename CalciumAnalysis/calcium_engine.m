%filename = 'day5_well1_control.tif';
% filename = '013 pintch 600g Z stack Time Series 1 day after 2nd injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-13.tif';
%filename = '012 pintch 300g Z stack Time Series 1 day after 2nd injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-12.tif';

%filename=['Images', filesep, '010 pintch 300g Z stack Time Series After 4th injection with Saline_1day_after_Pirt-GCaMP3 DRG live image-Orthogonal Projection-10.tif'];


%filename='011 100g stimulation Z stack Time Series 1 day after 4th injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-13.tif';
%filename='013 pintch 300g Z stack Time Series 1 day after 4th injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-14.tif';
%filename='014 pintch 600g Z stack Time Series 1 day after 4th injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-15.tif';
%filename='013 pintch 600g Z stack Time Series 1 day after 2nd injection with Saline_Pirt-GCaMP3 DRG live image-Orthogonal Projection-13.tif';
%filename='009 pintch 100g Z stack Time Series After 4th injection with Saline_1day_after_Pirt-GCaMP3 DRG live image-Orthogonal Projection-09.tif';
%filename='004 Large brush Burn 5D DRG L5 Spontaneous activity Time series Z stack Naive_Pirt-GCaMP3 DRG live image-Orthogonal Projection-04.tif';

% filename='006 small brush Burn 3D DRG L5 Spontaneous activity Time series Z stack Naive_Pirt-GCaMP3 DRG live image-Orthogonal Projection-06.tif';

function errorReport = calcium_engine(fileName, outputDirName)
errorReport = [];
dateCode='9-12-19';
timetotal = 15;
timestop = 15;
fprintf('Calling generateNPCGraph: %s\n', fileName);
[errorReport, cellinfo, framestop, bestsig1, bestsig2, badsig1, badsig2, adj, corrs, ...
    cutoff, ccov, bw, image, cellinfo_extended] = generateNPCgraph(fileName, timetotal, timestop);

if ~isempty(errorReport)
   return; 
end
%numCells = size(cellinfo, 1);
%scaleValues = false;
%generateDeltaPlot(cellinfo, 1:min(25, numCells), fileName, outputDirName, scaleValues);

% 2) Detect the size of the neurons according to the following criteria of sizes: <20um; >20 and <25; and >25;
% ZM : I believe that these have been captured at 2.5 um/pixel.
pixelwidth = 2.5; 
separations = [20 25]; % supply as um.
[cellSizeMetrics] = generateCellSizePlot(cellinfo_extended, pixelwidth, separations, fileName, outputDirName); 
disp(cellSizeMetrics)

% 
% 3) Count neurons showing calcium activity;
fprintf('%i of %i detected cells show calcium activity.\n', ... 
    size(cellinfo, 1), cellinfo_extended.totalCellCount);    

if size(cellinfo, 1) == 0
    fn = getFileName(fileName);
    errorReport = makeErrorStruct(['no active cells found in file ', fn], 1);
    return;
end

% generateCellClusters uses dbscan which was introduced in version R2019a.
% However, cytoNet uses an older version of Matlab
% 
% 4) Define ?cluster cells? (cells create a group of activation) according 
% to specific criteria like within 50um 3 or more cells activation;
% ZM : These are not great parameters for our method / test image combination, 
% so I used something else. This will definitely need to be tuned for each
% input and cell type. 
%idx = generateCellClusters(cellinfo, cellinfo_extended, 20, 4, pixelwidth, ...
%    image, fileName, outputDirName);  

% generateDeltaPlot uses stackedplot which was introduced in version R2018b.
% However, cytoNet uses an older version of Matlab
%
% 5) Report, for the neurons with calcium oscillation, the kinetic of the ?F/F0. 
%generateDeltaPlot(cellinfo, [10:10:100], fileName, outputDirName, true); 

generatePlots;

% Note that generatePlots changed the adjacency matrix, adj, so that low
% correlations are set to 0.  Now binarize adj.

adj = double(adj > 0);

% Mask is in bw

errorReport = writeGraphProperties(outputDirName, adj, bw, fileName);



% Adjacency Types
% 1: Centroid distance
% 2: Shared border
adjacencyType = 1;
threshold = 1.5;

if isempty(errorReport)
spatialGraph(outputDirName, bw, threshold, adjacencyType);
end

end
