function [GlobalMetrics, LocalMetrics, GlobalMetricNames, LocalMetricNames, processParams] = ...
    calculateNetworkMetrics(adjacencyMatrix)

tic;
processCompleted = true;

GlobalMetricNames = {'Node Count'; 'Edge Count'; 'Average Degree'; ... % basic network parameters
    'Network Density'; 'Variance in Degree'; 'Network Heterogeneity'; ... % degree-related metrics
    'Clustering Coefficient (normalized)'; 'Local Efficiency'; 'Average Neighbor Degree'; 'Variance in Neighbor Degree'; ...
    '4-star Motif Count'; '5-star Motif Count'; '6-star Motif Count'; 'Triangular Loop Count'; 'Pair Node Count'; 'Isolated Node Count'; ... % motif counts
    'Number of Connected Components'; 'Average Component Size (normalized)'; 'Variance in Component Size'; ... % metrics related to modularity
    'Network Diameter'; 'Network Efficiency (normalized)'; ... % geodesics
    'Assortativity'; 'Average Rich Club Metric'; 'Variance in Rich Club Metric'}; % degree-degree correlations

LocalMetricNames = {'Object Index'; ... % object identifiers
    'Degree'; 'Average Neighbor Degree'; 'Clustering Coefficient'; 'Local Efficiency'; ... % properties of local neighborhood
    'Closeness Centrality'; 'Betweenness Centrality'}; % centrality measures

% stop processing if number of nodes exceeds threshold
nNodes = numnodes(adjacencyMatrix);
if nNodes > 5000
    nEdges = -1;
    processCompleted = false;
    GlobalMetrics = 0;
    LocalMetrics = 0;
else
    
    cellIndex = 1:nNodes; % using MATLAB's indexing
    G = sparse(adjacencyMatrix);
    
    %% Degree-related metrics
        
    % Number of Edges
    nEdges = numedges(adjacencyMatrix);
    
    % stop processing if number of edges > 10000 or if total workspace
    % variable size exceeds 16GB
    if nEdges > 10000
        processCompleted = false;
        GlobalMetrics = 0;
        LocalMetrics = 0;
    elseif returnTotalVariableSize(whos) > 16
        processCompleted = false;
        fprintf('[calculateFunctionalNetworkMetrics] returnTotalVariableSize(whos)=%d\n', returnTotalVariableSize(whos));
        GlobalMetrics = 0;
        LocalMetrics = 0;
    else
        
        % Degree Distribution
        [kn, ~, ~] = degrees(adjacencyMatrix);
        kmax = max(kn); % maximum degree in network
        binCenters = 0:kmax;
        [degreeCounts, ~] = hist(kn, binCenters);
        
        % average degree (not normalized)
        avgeK = mean(kn);
        
        % variance of degree sequence normalized by maximum possible degree (n - 1)
        norm_kn = kn/(nNodes - 1);
        varK = var(norm_kn);
        
        % Average Neighborhood Degree for each Node
        NeighDegreeSequence = ave_neighbor_deg(adjacencyMatrix);
        
        % each element normalized by maximum possible neighborhood degree (n - 1)
        NormNeighDegreeSequence = NeighDegreeSequence/(nNodes - 1);
        
        % average of normalized neighbor degree sequence
        AvgeNeighborDegree = mean(NormNeighDegreeSequence);
        
        % variance of normalized neighbor degree sequence
        VarNeighborDegree = var(NormNeighDegreeSequence);
        
        % Network Density of graph = normalized average degree
        NetworkDensity = 2*nEdges/(nNodes*(nNodes - 1));
        
        % Network heterogeneity reflects tendency of hub nodes
        NetworkHeterogeneity = std(kn)/mean(kn);
        
        % Clustering Coefficient
        [C_local, C_local_all] = clust_coeff(adjacencyMatrix);
        C = mean(C_local); % global clustering coefficient
        
        % Local Efficiency
        [E_local_i, E_local_i_all] = local_efficiency(adjacencyMatrix);
        E_local = mean(E_local_i);
        
        % Assortativity
        Assortativity = pearson(adjacencyMatrix);
        
        % Rich-club metric for threshold degrees from 1 to maximum degree (n-1)
        RCM = zeros(1, (kmax-1));
        
        for i = 1:kmax-1
            RCM(i) = rich_club_metric(adjacencyMatrix, i);
        end
        
        %% Motif- and module-related metrics
        
        % number of nodes with no neighbors
        if numel(degreeCounts) >= 1
            %nSingleNodes = DegreeCounts(1)/nNodes;
            nSingleNodes = degreeCounts(1);
        else
            nSingleNodes = 0;
        end
        
        % number of independent pairs of nodes in graph
        if numel(degreeCounts) >= 2
            nPairNodes = 0.5*degreeCounts(2); %/nchoosek(nNodes, 2);
        else
            nPairNodes = 0;
        end
        
        % number of triangular loops
        if nNodes <=2
            nTriangularLoops = 0;
        else
            nTriangularLoops = trace(adjacencyMatrix^3)/6;
            %nTriangularLoops = nTriangularLoops/nchoosek(nNodes, 3);
        end
        
        % number of 4-node star motifs
        if nNodes <= 3
            nStar4 = 0;
        else
            nStar4 = num_star_motifs(adjacencyMatrix, 4); %/nchoosek(nNodes, 4);
        end
        
        % number of 5-node star motifs (normalized by n choose 5)
        if nNodes <= 4
            nStar5 = 0;
        else
            nStar5 = num_star_motifs(adjacencyMatrix, 5); %/nchoosek(nNodes, 5);
        end
        
        % number of 6-node star motifs (normalized by n choose 6)
        if nNodes <= 5
            nStar6 = 0;
        else
            nStar6 = num_star_motifs(adjacencyMatrix, 6); %/nchoosek(nNodes, 6);
        end
        
        % Number of connected components in the graph
        [ncc, ccLabels] = graphconncomp(G, 'directed', 'false');
        
        nccAbove1 = ncc; % number of connected components excluding single nodes
        ComponentSize = zeros(1, ncc);
        
        % removing single node and pair node connected components
        for c = 1:ncc
            ComponentSize(c) = numel(ccLabels(ccLabels == c));
            if ComponentSize(c) == 1
                nccAbove1 = nccAbove1 - 1;
                ccLabels(ccLabels == c) = 0; % all single node components labelled with 0
            end
        end
        
        nConnectedComponents = nccAbove1; % not normalized
        
        % finding component sizes
        UniqueComponentLabels = unique(ccLabels);
        nUniqueComponents = numel(UniqueComponentLabels);
        
        ComponentSize = zeros(1, nUniqueComponents-1);
        for i = 2:nUniqueComponents
            ccIndex = find(ccLabels == UniqueComponentLabels(i));
            % number of nodes and edges in connected component
            ComponentSize(i-1) = numel(ccIndex);
        end
        
        % size of components 2 nodes or larger, normalized by total nodes
        NormComponentSizeSequence = ComponentSize(ComponentSize >= 2)/nNodes;
        
        % average normalized component size
        AvgeComponentSizeNormalized = mean(NormComponentSizeSequence);
        
        % variance in normalized component size
        VarComponentSize = var(NormComponentSizeSequence);
        
        %% Computing centrality metrics
        
        ClosenessCentrality = closeness(G, adjacencyMatrix);
        
        if nNodes <= 2
            NodeBetweenness = zeros(1, nNodes);
            
        else
            NodeBetweenness = node_betweenness_faster(G, adjacencyMatrix);
        end
        
        %% Computing grodesic (path-length) metrics
        
        spl = graphallshortestpaths(G, 'directed', 'false');
        
        % Computing network diameter, normalized by longest path
        NetworkDiameter = max(spl(isfinite(spl)))/(nNodes-1);
        
        % Calculating path efficiency (inverse of path length)
        EfficiencyAllPaths = 1./spl;
        
        % efficiencies will be zeros only at diagonal
        E = mean(EfficiencyAllPaths(isfinite(EfficiencyAllPaths)));
        
        % NOTE: efficiency is not = 1/cpl because they are computed differently
        
        %% Generating 100 random graphs with the same number of nodes (N), edges (E), and same degree distribution (kn)
        
        nRandomGraphs = 100;
        
        C_random = zeros(nRandomGraphs, 1);
        GlobalEfficiencyRand_ri = zeros(nRandomGraphs, 1);
        E_local_random_ri = zeros(nRandomGraphs, 1);
        
        parfor ri = 1:nRandomGraphs
            AdjacencyMatrixRandom = random_graph(nNodes, [], [], 'sequence', kn);
            G_random = sparse(AdjacencyMatrixRandom);
            
            % Rich-club metric for threshold degrees from 1 to maximum degree (n-1)
            RCMr = zeros(1, (kmax-1));
            
            for i = 1:kmax-1
                RCMr(i) = rich_club_metric(AdjacencyMatrixRandom, i);
            end
            
            RCMrand_ri(ri, :) = RCMr;
            
            % calculating all shortest paths in random graph
            spl_random = graphallshortestpaths(G_random, 'directed', 'false');
            
            EfficiencyAllPaths = 1./spl_random;
            GlobalEfficiencyRand_ri(ri) = mean(EfficiencyAllPaths(isfinite(EfficiencyAllPaths))); % efficiencies will be zeros only at diagonal
            
            % Clustering Coefficient of random graph
            C_local_r = clust_coeff(AdjacencyMatrixRandom);
            
            C_random(ri) = mean(C_local_r);
            
            % Local efficiency of random graph
            E_local_r = local_efficiency(AdjacencyMatrixRandom);
            
            E_local_random_ri(ri) = mean(E_local_r);
            
        end
        
        % average rich-club metric sequence across 100 random graphs with minimum threshold degree = 1
        RCMrand = mean(RCMrand_ri);
        
        % normalized rich-club metric sequences
        RichClubMetricSequence = RCM./RCMrand;
        
        % average normalized rich-club metric across all degrees
        AvgeRichClubMetric = mean(RichClubMetricSequence);
        
        % variance in normalized rich-club metric sequence
        VarRichClubMetric = var(RichClubMetricSequence);
        
        % average global efficiency of 100 random graphs
        Erand = mean(GlobalEfficiencyRand_ri);
        
        % average clustering coefficient of 100 random graphs
        Crand = mean(C_random);
        
        % average local efficiency of 100 random graphs
        E_local_rand = mean(E_local_random_ri);
        
        Sigma = C/Crand; % normalized clustering coefficient
        
        EpsilonLocal = E_local/E_local_rand; % normalized local efficiency
        
        Epsilon = E/Erand; % normalized global efficiency
        
        %% Compiling network metric matrices
        
        LocalMetrics = [cellIndex', kn', NeighDegreeSequence', C_local_all', ...
            E_local_i_all', ClosenessCentrality, NodeBetweenness'];
        
        GlobalMetrics = [nNodes, nEdges, avgeK, ... % basic network parameters
            NetworkDensity, varK, NetworkHeterogeneity, ... % degree-related metrics
            Sigma, EpsilonLocal, AvgeNeighborDegree, VarNeighborDegree, ...
            nStar4, nStar5, nStar6, nTriangularLoops, nPairNodes, nSingleNodes, ... % motif counts
            nConnectedComponents, AvgeComponentSizeNormalized, VarComponentSize, ... % metrics related to modularity
            NetworkDiameter, Epsilon, ... % geodesics
            Assortativity, AvgeRichClubMetric, VarRichClubMetric]; % degree-degree correlations
        
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

