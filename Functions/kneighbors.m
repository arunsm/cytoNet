%% Function to compute the indices of k-neighborhood of a node 
%  (neighborhood of nodes that are 'k' links away from a given node)
%  http://strategic.mit.edu/docs/matlab_networks/kneighbors.m
function kneigh = kneighbors(adj, ind, k)
    adjk = adj;
    for i = 1:k-1 
        adjk = adjk.*adj; 
    end

    kneigh = find(adjk(ind, :) > 0); % returns indices of k-neighborhood
end