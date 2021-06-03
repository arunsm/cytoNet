% Computes the closeness centrality for every vertex: 1/sum(dist to all other nodes)
% For disconnected graphs can use: sum_over_t(2^-d(i,t)), idea Dangalchev (2006)
% C(i)=sum(2.^(-d)) if graph is disconnected, but sum w/o d(i)
% INPUTs: graph representation (adjacency matrix nxn)
% OUTPUTs: vector of centralities, nx1
% Source: social networks literature
% Other routines used: simple_dijkstra.m 
% GB, Last updated: October 9, 2009

function C=closeness(G, adj)

n = length(adj);
C=zeros(length(adj),1);  % initialize closeness vector

parfor i=1:length(adj); 
    %d = simple_dijkstra(adj,i); % distances from node i to every other node
    [d, ~, ~] = graphshortestpath(G, i, 'Directed' , 'false');
    d(i) = []; % d(i) is the distance from node i to itself
    C(i)=sum(1./2.^d); % using Dangalchev's version of centrality
end

C = C/(n-1); % normalizing