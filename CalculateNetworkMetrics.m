function [globalMetrics, localMetrics, globalMetricNames, localMetricNames, processParams] = calculateNetworkMetrics(adjacencyMatrix)

tic;
processCompleted = true;

globalMetricNames = {'Total Nodes'; 'Total Edges'; 'Average Degree'; ... % basic network parameters
    'Network Density'; 'Variance in Degree'; 'Network Heterogeneity'; ... % degree-related metrics
    'Clustering Coefficient (normalized)'; 'Local Efficiency'; 'Average Neighbor Degree'; 'Variance in Neighbor Degree'; ...
    '4-star Motif Count'; '5-star Motif Count'; '6-star Motif Count'; 'Triangular Loop Count'; 'Pair Node Count'; 'Isolated Node Count'; ... % motif counts
    'Number of Connected Components'; 'Average Component Size (normalized)'; 'Variance in Component Size'; ... % metrics related to modularity
    'Network Diameter'; 'Network Efficiency (normalized)'; ... % geodesics
    'Assortativity'; 'Average Rich Club Metric'; 'Variance in Rich Club Metric'}; % degree-degree correlations

localMetricNames = {'Object Index'; ... % object identifiers
    'Degree'; 'Average Neighbor Degree'; 'Clustering Coefficient'; 'Local Efficiency'; ... % properties of local neighborhood
    'Closeness Centrality'; 'Betweenness Centrality'}; % centrality measures

% stop processing if number of nodes exceeds threshold or if adjacency
% matrix is an empty matrix
nNodes = numnodes(adjacencyMatrix);
if nNodes > 5000 || nNodes <= 1
    nEdges = -1;
    globalMetrics = [];
    localMetrics = [];
    processCompleted = false;
    processParams.processCompleted = processCompleted;
    processParams.nNodes = nNodes;
    processParams.nEdges = nEdges;
    return;
else
    
    cellIndex = 1:nNodes; % using MATLAB's indexing
    G = sparse(adjacencyMatrix);
    
    %% Degree-related metrics
    
    nEdges = numedges(adjacencyMatrix); % number of edges
    
    % stop processing if number of edges > 10000 or if total workspace
    % variable size exceeds 16GB
    if nEdges > 10000
        processCompleted = false;
        processParams.processCompleted = processCompleted;
        processParams.nNodes = nNodes;
        processParams.nEdges = nEdges;
        globalMetrics = [];
        localMetrics = [];
        return;
    elseif returnTotalVariableSize(whos) > 16
        processCompleted = false;
        fprintf('[calculateNetworkMetrics] returnTotalVariableSize(whos)=%d\n', returnTotalVariableSize(whos));
        processParams.processCompleted = processCompleted;
        processParams.nNodes = nNodes;
        processParams.nEdges = nEdges;
        return;
    else
        
        % degree distribution
        [kn, ~, ~] = degrees(adjacencyMatrix); 
        kmax = max(kn);
        binCenters = 0:kmax;
        [degreeCounts, ~] = hist(kn, binCenters);
        
        % average degree
        avgeK = mean(kn);
        
        % variance of degree sequence normalized by maximum possible degree (n - 1)
        norm_kn = kn/(nNodes - 1);
        varK = var(norm_kn); 
        
        % average of normalized neighbor degree sequence
        neighDegreeSequence = ave_neighbor_deg(adjacencyMatrix);
        normNeighDegreeSequence = neighDegreeSequence/(nNodes - 1);
        avgeNeighborDegree = mean(normNeighDegreeSequence); 
        
        % variance of normalized neighbor degree sequence
        varNeighborDegree = var(normNeighDegreeSequence); 
        
        % network density of graph = normalized average degree
        networkDensity = 2*nEdges/(nNodes*(nNodes - 1)); 
        
        % network heterogeneity reflects tendency of hub nodes
        networkHeterogeneity = std(kn)/mean(kn);
        
        % clustering coefficient
        [C_local, C_local_all] = clust_coeff(adjacencyMatrix);
        clusteringCoefficient = mean(C_local);
        
        % local efficiency
        [E_local_i, E_local_i_all] = local_efficiency(adjacencyMatrix);
        localEfficiency = mean(E_local_i);
        
        % assortativity
        assortativity = pearson(adjacencyMatrix);
        
        % rich-club metric for threshold degrees from 1 to maximum degree (n-1)
        RCM = zeros(1, (kmax-1));
        for i = 1:kmax-1
            RCM(i) = rich_club_metric(adjacencyMatrix, i);
        end
        
        % average rich-club metric
        avgeRichClubMetric = mean(RCM);
        
        % variance in rich-club metric
        varRichClubMetric = var(RCM);
        
        %% motif- and module-related metrics
        
        % number of isolated nodes
        if numel(degreeCounts) >= 1
            nSingleNodes = degreeCounts(1);
        else
            nSingleNodes = 0;
        end
        
        % number of independent pairs of nodes in graph
        if numel(degreeCounts) >= 2
            nPairNodes = 0.5*degreeCounts(2);
        else
            nPairNodes = 0;
        end
        
        % number of triangular loops
        if nNodes <=2
            nTriangularLoops = 0;
        else
            nTriangularLoops = trace(adjacencyMatrix^3)/6;
        end
        
        % number of 4-node star motifs
        if nNodes <= 3
            nStar4 = 0;
        else
            nStar4 = num_star_motifs(adjacencyMatrix, 4);
        end
        
        % number of 5-node star motifs
        if nNodes <= 4
            nStar5 = 0;
        else
            nStar5 = num_star_motifs(adjacencyMatrix, 5);
        end
        
        % number of 6-node star motifs
        if nNodes <= 5
            nStar6 = 0;
        else
            nStar6 = num_star_motifs(adjacencyMatrix, 6);
        end
        
        % number of connected components excluding single nodes
        [ncc, ccLabels] = graphconncomp(G, 'directed', 'false');
        nccAbove1 = ncc; 
        componentSize = zeros(1, ncc);
        
        for c = 1:ncc
            componentSize(c) = numel(ccLabels(ccLabels == c));
            if componentSize(c) == 1
                nccAbove1 = nccAbove1 - 1;
                ccLabels(ccLabels == c) = 0;
            end
        end
        
        nConnectedComponents = nccAbove1;
        
        % average component size
        uniqueComponentLabels = unique(ccLabels);
        nUniqueComponents = numel(uniqueComponentLabels);
        
        componentSize = zeros(1, nUniqueComponents-1);
        for i = 2:nUniqueComponents
            ccIndex = find(ccLabels == uniqueComponentLabels(i));
            componentSize(i-1) = numel(ccIndex);
        end
        
        componentSizeSequence = componentSize(componentSize >= 2);
        avgeComponentSize = mean(componentSizeSequence);
        
        % variance in component size
        varComponentSize = var(componentSizeSequence);
        
        %% centrality metrics
        
        closenessCentrality = closeness(G, adjacencyMatrix);
        
        if nNodes <= 2
            nodeBetweenness = zeros(1, nNodes); 
        else
            nodeBetweenness = node_betweenness_faster(G, adjacencyMatrix);
        end
        
        %% geodesic metrics
        
        spl = graphallshortestpaths(G, 'directed', 'false');
        
        % network diameter
        networkDiameter = max(spl(isfinite(spl)));
        
        % path efficiency (inverse of path length)
        efficiencyAllPaths = 1./spl;
        globalEfficiency = mean(efficiencyAllPaths(isfinite(efficiencyAllPaths)));
       
        %% compile network metric matrices
        
        localMetrics = [cellIndex', kn', neighDegreeSequence', C_local_all', ...
            E_local_i_all', closenessCentrality, nodeBetweenness'];
        
        globalMetrics = [nNodes, nEdges, avgeK, ... % basic network parameters
            networkDensity, varK, networkHeterogeneity, ... % degree-related metrics
            clusteringCoefficient, localEfficiency, avgeNeighborDegree, varNeighborDegree, ...
            nStar4, nStar5, nStar6, nTriangularLoops, nPairNodes, nSingleNodes, ... % motif counts
            nConnectedComponents, avgeComponentSize, varComponentSize, ... % metrics related to modularity
            networkDiameter, globalEfficiency, ... % geodesics
            assortativity, avgeRichClubMetric, varRichClubMetric]; % degree-degree correlations
        
    end
end

totalWorkspaceSize = returnTotalVariableSize(whos);
processingTime = toc;

processParams.processCompleted = processCompleted;
processParams.totalWorkspaceSize = totalWorkspaceSize;
processParams.processingTime = processingTime;
processParams.nNodes = nNodes;
processParams.nEdges = nEdges;
end

