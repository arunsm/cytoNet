function Mask2GraphStruct = createMask2GraphStruct(Mask2GraphRequested, adjacencyType, threshold, maskLayerNumber)
   
    % boolean value indicating whether graph implementation is required
    Mask2GraphStruct.Mask2GraphRequested = Mask2GraphRequested; 

    % numeric indicator for type of threshold used to build graph
    % 1 = shared border pixels
    % 2 = centroid-centroid distance
    Mask2GraphStruct.adjacencyType = adjacencyType;
    
    % threshold value for determining adjacency between objects
    Mask2GraphStruct.threshold = threshold;
    
    % layer number to be used for creating graph
    Mask2GraphStruct.maskLayerNumber = maskLayerNumber;

end