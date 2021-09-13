
function [errorReport, imageProcessingParams] = Mask2Graph(maskInputPath, outputPath, adjacencyType, threshold, maskLayerNumber, errorReport)

s = createImageProcessingParams(0, 0, 0, 0, 0, 0, 0, 0);
imageProcessingParams = repmat(s, 0);

flag = false;

%% starting parallel pool
myCluster = parcluster('local');
myCluster.JobStorageLocation = pwd;
poolObj = parpool(myCluster);

%% Loading parameters and initializing variables

% read all files in folder specified by maskInputPath
maskReadPath = strcat(maskInputPath, filesep, '*');
d = dir(maskReadPath);
fnme = {d.name};

% accept only files and not folders
maskIdx = cellfun(@isdir, fnme);
fnmeMasks = fnme(~maskIdx);
nMaskFiles = numel(fnmeMasks);

%% process all images in folder
for i = 1:nMaskFiles
    writeLog(sprintf('[Mask2Graph] Mask number %d/%d', i, nMaskFiles));
    
    currentMaskFile = strcat(maskInputPath, filesep, fnmeMasks{i});
    currentMaskBaseName = createBaseName(fnmeMasks{i});
    
    % finding layer number and base file name for already segmented images
    idxLayer = strfind(currentMaskBaseName, 'layer');
    idxMASKS = strfind(currentMaskBaseName, 'MASKS');
    
    s = '';
    if isempty(idxLayer)
        layerNumberParam = '';
    else
        layerNumberParam = currentMaskBaseName(idxLayer+5);
        s = strcat('-layer', layerNumberParam);
        
        if ~isempty(idxMASKS)
            s = strcat(s, '-MASKS');
        end
    end
    
    originalImageName = strrep(currentMaskBaseName, s, '');
    
    % make sure currentMaskFile is a readable image file
    try
        info = imfinfo(currentMaskFile);
    catch
        errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is not a readable image file'], 1);
        continue;
    end
    
    imageDimensions = [info.Width info.Height];
    nLayers = numel(info); % number of layers in image
    
    % maskLayerNumber will be an empty vector if segmentation has been requested
    if isempty(maskLayerNumber)
        maskLayerNumber = 1;
    end
    
    for k = 1:numel(maskLayerNumber) % loop over number of layers in image to be processed
        writeLog(sprintf('[Mask2Graph] Mask layer %d/%d', k, numel(maskLayerNumber)));
        if maskLayerNumber == 1
            writeLog('[Mask2Graph] maskLayerNumber is 1');
            currentMask = imread(currentMaskFile);
            layerSuffix = '';
        else
            writeLog('[Mask2Graph] maskLayerNumber is NOT 1');
            if nLayers<maskLayerNumber(k)
                writeLog('[Mask2Graph] nlayers<maskLayerNumber !!!');
                errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' does not have ', num2str(maskLayerNumber(k)), ' layers'], 1);
                continue;
            end
            
            currentMask = imread(currentMaskFile, maskLayerNumber(k));
            layerSuffix = strcat('layer', num2str(maskLayerNumber(k)));
            layerNumberParam = num2str(maskLayerNumber(k));
        end
        writeLog('[Mask2Graph] Current mask read');
        
        flag = true; % this means we have usable input
        
        % if image is RGB, convert to grayscale
        if size(currentMask, 3) > 1
            currentMask = rgb2gray(currentMask(:, :, 1:3));
            errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is an RGB image: converting to grayscale'], 1);
        end
        writeLog('[Mask2Graph] Converted to grayscale if necessary');
        
        % make sure image is logical
        if ~islogical(currentMask)
            currentMask = currentMask > 0;
            errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is not a binary image: converting to binary'], 1);
        end
        
        writeLog('[Mask2Graph] Checked for logical values');
        writeLog('[Mask2Graph] fprinted');
        
        % call generateSpatialGraph to generate spatial graph
        cellInfoAllCells = generateSpatialGraph(currentMask, threshold, adjacencyType);
        writeLog('[Mask2Graph] Completed generateSpatialGraph');
        
        currentMaskName = strcat(currentMaskBaseName, layerSuffix);
        writeLog('[Mask2Graph] Set currentMaskName');
        cellInfoAllCells.fileName = currentMaskName;
        
        % call calculateNetworkMetrics to compute graph-based metrics
        [cellInfoAllCells.globalMetrics, cellInfoAllCells.localMetrics, ...
            cellInfoAllCells.globalMetricNames, cellInfoAllCells.localMetricNames, processParams] = ...
            calculateNetworkMetrics(cellInfoAllCells.adjacencyMatrixBinary);
        writeLog('[Mask2Graph] Completed CalculateNetworkMetrics');
        
        % calculate all network metrics for null model
        cellInfoAllCells.globalMetricsRandom = calculateNetworkMetricsRandom(cellInfoAllCells);
        
        writeLog(sprintf('[Mask2Graph] WriteNetworkMetricsFlag=%d', flag));
        % plot processed image and write network metrics to file
        if processParams.processCompleted
            plotProcessedImage(outputPath, cellInfoAllCells);
            writeNetworkMetrics(outputPath, cellInfoAllCells);
        else
            errorReport(end+1) = makeErrorStruct(['processing for file ', fnmeMasks{i}, ' exceeds computing capacity; graph not created'], 1);
        end
        
        imageProcessingParams(end+1) = createImageProcessingParams(originalImageName, ...
            layerNumberParam, imageDimensions, processParams.nNodes, processParams.nEdges, ...
            processParams.totalWorkspaceSize, processParams.processingTime, processParams.processCompleted);
    end
end

%% clean up
writeLog('[Mask2Graph] Deleting pool object');
delete(poolObj);
writeLog('[Mask2Graph] Done');

end
