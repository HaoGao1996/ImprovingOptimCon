function [weightednet, downnodes] = get_weightednet(A, gene_exp, der_gene, epsilon)
% FUNCTION: get the weighted network of downstream nodes
% Input:
% A is the adjacency matrix;
% gene_exp is the expression data of genes in GRN;
% der_gene is the deregulationg gene;
% epsilon is the value correpondings to a p-value of 0.05 based on the
% theoretical distribution of the PCC.
% Output:
% weightednet is the weighted adjacency matrix of downstream nodes;
% downnodes is the set of downstream nodes of deregulation genes; 

%% get downrelations
downrelation=[];
downweight=[];
downnodes = der_gene;
Q = der_gene; % quene for BFS
while ~isempty(Q)
    x = Q(1);
    
    [relation, weight, x_downnode] = get_downrelation(A, x, der_gene, gene_exp, epsilon);
    downrelation = [downrelation; relation]; 
    downweight = [downweight; weight];
    
    x_downnode(find(ismember(x_downnode, downnodes))) = []; % delete exist nodes
    
    downnodes = [downnodes; x_downnode];
    Q = [Q; x_downnode]; % put new nodes into Quene
    
    Q(1) = []; % delete the head node
end

%% construct weighted network
num = length(downnodes);
weightednet = zeros(num);
[~, loc_idx] = ismember(downrelation, downnodes);
for i = 1:size(downrelation, 1)
    weightednet(loc_idx(i, 1), loc_idx(i, 2)) = downweight(i); % A12 <=> 1->2
end

end