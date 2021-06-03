%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION:  Function to compute image metrics based on
%               network analysis
%
% OUTPUT:       Global metrics, Local Metrics, Adjacency Matrix and Cell
%               Locations
%
% AUTHOR:       Arun Mahadevan
% LAST UPDATED: 6/18/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GlobalMetrics, LocalMetrics, AdjacencyMatrix, MaskBoundaries, ...
    CellLocations, GlobalMetricNames, LocalMetricNames, processParams] = ...
    CalculateNetworkMetrics(Masks, threshold, adjacencyType)

tic;
processCompleted = true;

GlobalMetricNames = {'Node Count'; 'Edge Count'; 'Fraction Area Cells'; 'Average Degree'; ... % basic network parameters
    'Network Density'; 'Variance in Degree'; 'Network Heterogeneity'; ... % degree-related metrics
    'Clustering Coefficient (normalized)'; 'Local Efficiency'; 'Average Neighbor Degree'; 'Variance in Neighbor Degree'; ...
    '4-star Motif Count'; '5-star Motif Count'; '6-star Motif Count'; 'Triangular Loop Count'; 'Pair Node Count'; 'Isolated Node Count'; ... % motif counts
    'Number of Connected Components'; 'Average Component Size (normalized)'; 'Variance in Component Size'; ... % metrics related to modularity
    'Network Diameter'; 'Network Efficiency (normalized)'; ... % geodesics
    'Assortativity'; 'Average Rich Club Metric'; 'Variance in Rich Club Metric'}; % degree-degree correlations

LocalMetricNames = {'Object Index'; 'Centroid(x)'; 'Centroid(y)'; ... % object identifiers
    'Degree'; 'Average Neighbor Degree'; 'Clustering Coefficient'; 'Local Efficiency'; ... % properties of local neighborhood
    'Closeness Centrality'; 'Betweenness Centrality'}; % centrality measures

% creating a connected component structure to represent masks
CC = bwconncomp(Masks);
nNodes = CC.NumObjects;
ImageDimensions = CC.ImageSize;

% finding object boundaries
MaskBoundaries = bwboundaries(Masks, 'noholes');

if nNodes > 5000
    
    nEdges = -1;
    processCompleted = false;
    GlobalMetrics = 0;
    LocalMetrics = 0;
    AdjacencyMatrix = 0;
    CellLocations = [0 0];
    
else
    
    % Extracting cell centroids
    c = regionprops(CC, 'Centroid');
    CellLocations = cat(1, c.Centroid);
    
    % Computing equivalent diameters of cells
    d = regionprops(CC, 'EquivDiameter');
    diameters = cat(1, d.EquivDiameter);
    
    CellIndex = 1:nNodes; % using MATLAB's indexing
    
    %% Creating network representation using centroid-centroid proximity
    
    if adjacencyType == 1
        
        % Setting distance threshold for each pair of nodes as the average
        % diameter of that pair of cells multiplied by scaling factor
        % (threshold)
        distanceThreshold = threshold*squareform(pdist(diameters, @averageDiameter));
        
        %distanceThreshold = 282; % hard threshold option?
        
        % Computing adjacency matrix
        weights = squareform(pdist(CellLocations, 'Euclidean'));
        AdjacencyMatrix = double(weights < distanceThreshold);
        
    end
    
    %% Adjacency matrix based on shared border pixels
    
    if adjacencyType == 2
        
        % Splitting each cell mask into a separate image in a cell array
        MasksCellArray = cell(nNodes, 1);
        for i = 1:nNodes
            mask_i = false(ImageDimensions);
            mask_i(CC.PixelIdxList{i}) = true;
            MasksCellArray{i} = mask_i;
        end
        
        allMasksExtended = cell([nNodes 1]); % initialize new binary image containing all extended masks
        filter = [0 1 0; 1 1 1; 0 1 0]; % filter for mask extension
        
        % extending masks out by 2 pixels
        for i = 1:nNodes
            mask_i = MasksCellArray{i};
            mask_i_extended_1pixel = logical(imfilter(mask_i, filter)); % extending mask i based on filter
            mask_i_extended_2pixels = logical(imfilter(mask_i_extended_1pixel, filter));
            
            allMasksExtended{i} = mask_i_extended_2pixels;
        end
        
        BorderLength = squareform(sharedPixels(allMasksExtended));
        
        AdjacencyMatrix = double(BorderLength > threshold);
        
    end
    
    %% Assigning diagonal elements of AdjacencyMatrix as zero
    for i = 1:nNodes
        AdjacencyMatrix(i, i) = 0;
    end
    
    G = sparse(AdjacencyMatrix);
    
    %% Computing degree-related metrics
    
    % Percent cell area covered
    TotalAreaCells = 0;
    A = regionprops(CC, 'Area');
    for i = 1:nNodes
        TotalAreaCells = TotalAreaCells + A(i).Area;
    end
    PercentAreaCells = TotalAreaCells/(ImageDimensions(1)*ImageDimensions(2));
    
    % Number of Edges
    nEdges = sum(AdjacencyMatrix(:))/2;
    
    % stopping processing if number of edges > 10000
    if nEdges > 10000
        processCompleted = false;
        GlobalMetrics = 0;
        LocalMetrics = 0;
        
    % stopping processing if total workspace variable size exceeds 16GB
    elseif returnTotalVariableSize(whos) > 16
        processCompleted = false;
        GlobalMetrics = 0;
        LocalMetrics = 0;
        
    else
        
        % Degree Distribution
        [kn, ~, ~] = degrees(AdjacencyMatrix);
        
        % Computing Degree Distribution
        kmax = max(kn); % Maximum degree in network
        binCenters = 0:kmax;
        [DegreeCounts, ~] = hist(kn, binCenters);
        
        % average degree (not normalized)
        avgeK = mean(kn);
        
        % variance of degree sequence normalized by maximum possible degree (n - 1)
        norm_kn = kn/(nNodes - 1);
        varK = var(norm_kn);
        
        % Average Neighborhood Degree for each Node
        NeighDegreeSequence = ave_neighbor_deg(AdjacencyMatrix);
        
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
        [C_local, C_local_all] = clust_coeff(AdjacencyMatrix);
        C = mean(C_local); % global clustering coefficient
        
        % Local Efficiency
        [E_local_i, E_local_i_all] = local_efficiency(AdjacencyMatrix);
        E_local = mean(E_local_i);
        
        % Assortativity
        Assortativity = pearson(AdjacencyMatrix);
        
        % Rich-club metric for threshold degrees from 1 to maximum degree (n-1)
        RCM = zeros(1, (kmax-1));
        
        for i = 1:kmax-1
            RCM(i) = rich_club_metric(AdjacencyMatrix, i);
        end
        
        %% Computing motif- and module-related metrics
        
        % number of nodes with no neighbors
        if numel(DegreeCounts) >= 1
            %nSingleNodes = DegreeCounts(1)/nNodes;
            nSingleNodes = DegreeCounts(1);
        else
            nSingleNodes = 0;
        end
        
        % number of independent pairs of nodes in graph
        if numel(DegreeCounts) >= 2
            nPairNodes = 0.5*DegreeCounts(2); %/nchoosek(nNodes, 2);
        else
            nPairNodes = 0;
        end
        
        % number of triangular loops
        if nNodes <=2
            nTriangularLoops = 0;
        else
            nTriangularLoops = trace(AdjacencyMatrix^3)/6;
            %nTriangularLoops = nTriangularLoops/nchoosek(nNodes, 3);
        end
        
        % number of 4-node star motifs
        if nNodes <= 3
            nStar4 = 0;
        else
            nStar4 = num_star_motifs(AdjacencyMatrix, 4); %/nchoosek(nNodes, 4);
        end
        
        % number of 5-node star motifs (normalized by n choose 5)
        if nNodes <= 4
            nStar5 = 0;
        else
            nStar5 = num_star_motifs(AdjacencyMatrix, 5); %/nchoosek(nNodes, 5);
        end
        
        % number of 6-node star motifs (normalized by n choose 6)
        if nNodes <= 5
            nStar6 = 0;
        else
            nStar6 = num_star_motifs(AdjacencyMatrix, 6); %/nchoosek(nNodes, 6);
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
        
        ClosenessCentrality = closeness(G, AdjacencyMatrix);
        
        if nNodes <= 2
            NodeBetweenness = zeros(1, nNodes);
            
        else
            NodeBetweenness = node_betweenness_faster(G, AdjacencyMatrix);
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
        
        LocalMetrics = [CellIndex', CellLocations, kn', NeighDegreeSequence', C_local_all', ...
            E_local_i_all', ClosenessCentrality, NodeBetweenness'];
        
        GlobalMetrics = [nNodes, nEdges, PercentAreaCells, avgeK, ... % basic network parameters
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

