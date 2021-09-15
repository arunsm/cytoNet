function globalMetricsRandom = calculateNetworkMetricsRandom(cellInfoAllCells)

rng('default');
nRandomGraphs = 100;

globalMetricNames = {'Total Nodes'; 'Total Edges'; 'Average Degree'; ... % basic network parameters
    'Network Density'; 'Variance in Degree'; 'Network Heterogeneity'; ... % degree-related metrics
    'Clustering Coefficient (normalized)'; 'Local Efficiency'; 'Average Neighbor Degree'; 'Variance in Neighbor Degree'; ...
    '4-star Motif Count'; '5-star Motif Count'; '6-star Motif Count'; 'Triangular Loop Count'; 'Pair Node Count'; 'Isolated Node Count'; ... % motif counts
    'Number of Connected Components'; 'Average Component Size (normalized)'; 'Variance in Component Size'; ... % metrics related to modularity
    'Network Diameter'; 'Network Efficiency (normalized)'; ... % geodesics
    'Assortativity'; 'Average Rich Club Metric'; 'Variance in Rich Club Metric'}; % degree-degree correlations

nGlobalMetrics = numel(globalMetricNames);
globalMetricsRandom = zeros(nRandomGraphs, nGlobalMetrics);

adjacencyType = cellInfoAllCells.adjacencyType;
if adjacencyType == 2 % if spatial graph, use spatial null model (only type 2 graphs)
    mask = cellInfoAllCells.mask;
    threshold = cellInfoAllCells.threshold;
    parfor i = 1:nRandomGraphs
        [~, adjacencyMatrixRandom, ~] = spatialNullModel(mask, threshold, adjacencyType);
        [globalMetricsRandom(i, :), ~, ~, ~, ~] = calculateNetworkMetrics(adjacencyMatrixRandom);
    end
    
else % if functional graph or type 1 spatial graph, use random rewiring null model
    nNodes = cellInfoAllCells.nNodes;
    [kn, ~, ~] = degrees(cellInfoAllCells.adjacencyMatrixBinary);
    parfor i = 1:nRandomGraphs
        adjacencyMatrixRandom = random_graph(nNodes, [], [], 'sequence', kn);
        [globalMetricsRandom(i, :), ~, ~, ~, ~] = calculateNetworkMetrics(adjacencyMatrixRandom)
    end
end

end