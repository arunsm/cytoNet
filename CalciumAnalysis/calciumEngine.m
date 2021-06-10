function errorReport = calciumEngine(fileName, outputDirName)

errorReport = [];

[errorReport, cellInfo, cellInfoAllCells, mask, maxImage, functionalAdjacencyMatrix] = generateFunctionalGraph(filename);

if ~isempty(errorReport)
    return;
end

if size(cellInfo, 1) == 0
    fn = getFileName(fileName);
    errorReport = makeErrorStruct(['no active cells found in file ', fn], 1);
    return;
end

generatePlots;

errorReport = writeGraphProperties(outputDirName, functionalAdjacencyMatrix, mask, fileName);

% create spatial graph
adjacencyType = 1;
threshold = 1.5;

if isempty(errorReport)
    spatialGraph(outputDirName, bw, threshold, adjacencyType);
end

end