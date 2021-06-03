function ImageProcessingParamsStruct = createImageProcessingParams(imageName, ...
    layerNumber, imageDimensions, nNodes, nEdges, totalWorkspaceSize, processingTime, processCompleted)

    % Image base name
    ImageProcessingParamsStruct.imageName = imageName; 
    
    % Layer number for multi-layered images
    ImageProcessingParamsStruct.layerNumber = layerNumber; 
    
    % Image dimensions as widthxheight in pixels
    ImageProcessingParamsStruct.imageDimensions = imageDimensions; 

    % Number of nodes counted from binary image
    ImageProcessingParamsStruct.nNodes = nNodes;
    
    % Number of edges in graph
    ImageProcessingParamsStruct.nEdges = nEdges;
    
    % Total workspace size used by CalculateNetworkMetrics.m (in gigabytes)
    ImageProcessingParamsStruct.totalWorkspaceSize = totalWorkspaceSize;

    % Total processing time used by CalculateNetworkMetrics.m (in seconds)
    ImageProcessingParamsStruct.processingTime = processingTime;
    
    % Boolean indicator for successful completion of graph analysis
    ImageProcessingParamsStruct.processCompleted = processCompleted;
    
end