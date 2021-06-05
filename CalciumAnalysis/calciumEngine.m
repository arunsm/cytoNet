function errorReport = calciumEngine(fileName, outputDirName)

errorReport = [];
timetotal = 15;
timestop = 15;

[errorReport, cellinfo, framestop, bestsig1, bestsig2, badsig1, badsig2, adj, corrs, ...
    cutoff, ccov, bw, image, cellinfo_extended] = generateNPCgraph(fileName, timetotal, timestop);

if ~isempty(errorReport)
    return;
end

if size(cellinfo, 1) == 0
    fn = getFileName(fileName);
    errorReport = makeErrorStruct(['no active cells found in file ', fn], 1);
    return;
end

generatePlots;

errorReport = writeGraphProperties(outputDirName, adj, bw, fileName);

% create spatial graph
adjacencyType = 1;
threshold = 1.5;

if isempty(errorReport)
    spatialGraph(outputDirName, bw, threshold, adjacencyType);
end

end