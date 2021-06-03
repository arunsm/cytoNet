% Computes clustering coefficient of each node in graph
% INPUT: adjacency matrix representation of a graph
% OUTPUT: graph average clustering coefficient and clustering coefficient
% Other routines used: degrees.m, kneighbors.m, subgraph.m

function [C_truncated, C] = clust_coeff(adj)

n = length(adj);
[deg, ~, ~] = degrees(adj);
C = zeros(1, n); % initialize clustering coefficient
ind = ones(1, n);

for i = 1:n
  
  if deg(i) <= 1 
      ind(i) = 0;
      continue; 
  end 
  
  neigh=kneighbors(adj,i,1);
  sub_neigh = subgraph(adj, neigh); % subgraph of neighborhood
  edges_s = sum(sub_neigh(:))/2; % number of edges in neighborhood
   
  C(i)=2*edges_s/deg(i)/(deg(i)-1);

end

C_truncated = C(ind == 1); % only clustering coefficients of nodes with degree > 1