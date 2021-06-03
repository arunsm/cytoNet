function Mask2GraphStruct = createMask2GraphStruct(Mask2GraphRequested, thresholdType, threshold, maskLayerNumber)
   
    % boolean value indicating whether graph implementation is required
    Mask2GraphStruct.Mask2GraphRequested = Mask2GraphRequested; 

    % numeric indicator for type of threshold used to build graph
    % 1 = centroid-centroid distance
    % 2 = shared border pixels
    Mask2GraphStruct.thresholdType = thresholdType;
    
    % threshold value for determining adjacency between objects
    Mask2GraphStruct.threshold = threshold;
    
    % layer number to be used for creating graph
    Mask2GraphStruct.maskLayerNumber = maskLayerNumber;

end