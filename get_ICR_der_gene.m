function ICR_der_gene = get_ICR_der_gene(weightednet, downnodes, der_gene, lambda)
% FUNCTION: get the indirect control region of a deregulation gene
% Input:
% weightednet is the weighted adjacency matrix of downstream nodes;
% downnodes is the set of downstream nodes of deregulation genes; 
% lambda is the threshold of ICV
% Output: 
% ICR_i: indirect control region of a deregulation gene i

weightednet = sparse(weightednet);
num = length(weightednet); % length of the weightednet
u = find(downnodes == der_gene); % the lable of der_gene in weightednet

% d shortest distance
% pred is the shortest path
[d, pred] = dijkstra_sp(weightednet, u);

ICV_der_gene = 1./d;
ICR_der_gene = [];
for i = 1:num
    if ICV_der_gene(i) >= lambda % 如果ICV的值大于lambda = 0.3就认为是indirect control region
        SP = get_shortest_path(u, i, pred);
        ICR_der_gene = [ICR_der_gene; downnodes(SP)];
    end
end

ICR_der_gene = unique(ICR_der_gene);
end


