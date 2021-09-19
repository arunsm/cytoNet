function segmentationStruct = createSegmentationStruct(segmentationRequested, segmentationType, segmentationLayerNumber)
    segmentationStruct = struct;
    
    % boolean value indicating whether segmentation is required
    segmentationStruct.segmentationRequested = segmentationRequested; 
    
    % string indicating type of segmentation; current options are
    % 't' - locally adaptive thresholding
    % 'tw' - locally adaptive thresholding followed by watershed
    segmentationStruct.segmentationType = segmentationType; 
    
    % layer number to be used for segmentation (for layered .tif files)
    segmentationStruct.segmentationLayerNumber = segmentationLayerNumber;
end