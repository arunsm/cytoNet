function [errorReport, imageProcessingParams] = Engine(InputPath, OutputPath, segmentationStruct, Mask2GraphStruct)

s = makeErrorStruct('', '');
errorReport = repmat(s, 0);

segmentationRequested = segmentationStruct.segmentationRequested;
segmentationType = segmentationStruct.segmentationType;
segmentationLayerNumber = segmentationStruct.segmentationLayerNumber;

Mask2GraphRequested = Mask2GraphStruct.Mask2GraphRequested;
adjacencyType = Mask2GraphStruct.adjacencyType;
threshold = Mask2GraphStruct.threshold;
maskLayerNumber = Mask2GraphStruct.maskLayerNumber;

if segmentationRequested
    
    % make new folder to store masks
    maskInputPath = 'MASKS';
    mkdir(maskInputPath);
    
    % read all files in folder specified by InputPath
    ImageReadPath = strcat(InputPath, filesep, '*');
    d = dir(ImageReadPath);
    fnme = {d.name};
    
    % accept only files and not folders
    imageIdx = cellfun(@isdir, fnme);
    fnmeImages = fnme(~imageIdx);
    
    nImages = numel(fnmeImages);
    singleCellMetrics = {}; % cell array containing cell shape and intensity metrics
    singleCellMetricNames = {}; % cell array containing single-cell metric names
    ctr = 1;
    for i = 1:nImages
        currentImageFile = strcat(InputPath, filesep, fnmeImages{i});
        
        writeLog(sprintf('[Engine] Segmenting %s', currentImageFile));
        
        % make sure currentImageFile is a readable image file
        try
            info = imfinfo(currentImageFile);
        catch
            errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is not a readable image file; skipping this file'], 1);
            continue;
        end
        
        % using segmentationLayerNumber to extract relevant images for
        % segmentation
        nLayers = numel(info); % total number of layers in image
        for k = 1:numel(segmentationLayerNumber)
            
            writeLog(sprintf('[Engine] Getting image number %d', k));
            
            if nLayers == 1 || all(segmentationLayerNumber == 1) % nLayers will be 1 if image is a non-layered image
                currentImage = imread(currentImageFile);
                layerSuffix = '';
                
                if nLayers<segmentationLayerNumber(k)
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' does not have ', num2str(segmentationLayerNumber(k)), ' layers; skipping this file'], 1);
                    continue;
                end
            else
                if nLayers<segmentationLayerNumber(k)
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' does not have ', num2str(segmentationLayerNumber(k)), ' layers; skipping this file'], 1);
                    continue;
                end
                
                currentImage = imread(currentImageFile, segmentationLayerNumber(k));
                layerSuffix = strcat('-layer', num2str(k));
            end
            
            % if image is RGB, converting to grayscale
            if size(currentImage, 3) >= 3
                currentImage = rgb2gray(currentImage(:, :, 1:3));
                errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is an RGB image: converting to grayscale'], 1);
            else
                % Check for unusual situation of a two channel image
                if size(currentImage, 3) == 2
                    currentImage = currentImage(:, :, 1);
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' has only two channels; using only the first channel'], 1);
                end
            end
            
            currentImageBaseName = createBaseName(fnmeImages{i});
            currentImageName = strcat(currentImageBaseName, layerSuffix);
            
            % performing segmentation if image is non-logical
            if ~islogical(currentImage)
                
                fprintf('\tperforming segmentation for image %s\n', fnmeImages{i});
                
                switch segmentationType
                    case 't'
                        Masks = segmentThreshold(currentImage);
                    case 'tw'
                        Masks = segmentWatershed(currentImage);
                    otherwise
                        error(['[Engine] Unexpected segmentation type: ', segmentationType]);
                end
                
                singleCellMetrics{ctr, 1} = currentImageName;
                
                % calculating single-cell shape and intensity metrics
                [singleCellMetricNames, singleCellMetrics{ctr, 2}] = calculateSingleCellMetrics(currentImage, Masks);
                
            else
                errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is a binary image: did not perform segmentation'], 1);
                Masks = currentImage;
            end
            
            % writing mask files to mask folder
            
            writeLog('[Engine] Writing mask file');
            
            MaskFileName = strcat(maskInputPath, filesep, currentImageBaseName, layerSuffix, '-MASKS.tif');
            imwrite(Masks, MaskFileName);
            
            ctr = ctr + 1;
        end
    end
    
    % writing single-cell metrics to .csv files
    WriteSingleCellMetrics(OutputPath, singleCellMetricNames, singleCellMetrics);
    
else
    maskInputPath = InputPath;
end

if Mask2GraphRequested
    
    writeLog('[Engine] Calling Mask2Graph');
    
    
    % call Mask2Graph to convert mask info to graph output
    [errorReport, imageProcessingParams] = Mask2Graph(maskInputPath, OutputPath, adjacencyType, threshold, maskLayerNumber, errorReport);
    
else
    
    imageProcessingParams = struct; % empty struct indicates no graph analysis was done
    
    % transfer masks to output directory
    maskFiles = strcat(maskInputPath, filesep, '*-MASKS.tif');
    d = dir(maskFiles);
    
    if numel(d) > 0
        movefile(maskFiles, OutputPath);
    end
end



writeLog('[Engine] Done');


