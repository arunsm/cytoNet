function [LE_truncated, LE] = local_efficiency(adj)

n = length(adj);
LE = zeros(1, n);
ind = ones(1, n);
[deg, ~, ~] = degrees(adj);

for i = 1:n
  
  neigh=kneighbors(adj,i,1);
  
  if(deg(i)<=1)
      ind(i) = 0;
      continue;
  end
  
  sub_neigh = subgraph(adj, neigh); % subgraph of neighborhood
   
  G_sub = sparse(sub_neigh);
  
  spl = graphallshortestpaths(G_sub, 'directed', 'false');
  
  EfficiencyAllPaths = 1./spl;
    
  LE(i) = mean(EfficiencyAllPaths(isfinite(EfficiencyAllPaths)));

end

LE_truncated = LE(ind == 1); % only considering nodes with degree > 1