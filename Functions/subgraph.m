%% Function to return sub-graph given specific node indices
%  http://strategic.mit.edu/docs/matlab_networks/subgraph.m
function adj_sub = subgraph(adj, S)
    adj_sub = adj(S,S);
end