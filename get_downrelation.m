function [relation, weight, x_downnode] = get_downrelation(A, x, der_gene, gene_exp, epsilon)
% FUNCTION: get the weighted network of downstream nodes
% Input:
% A is the adjacency matrix;
% x is the source node;
% gene_exp is the expression data of genes in GRN;
% der_gene is the deregulationg gene;
% epsilon is the value correpondings to a p-value of 0.05 based on the
% theoretical distribution of the PCC.
% Output:
% relation is set of edge from x;
% weight is weight of the set of relations; 
% x_downnode is the set of nodes pointed by node x;

x_downnode = find(A(:,x));
relation = [];
weight = [];

if ~isempty(x_downnode)
    for i = 1:length(x_downnode)
        relation = [relation; [x, x_downnode(i)]];
        w = comp_w(x, x_downnode(i),der_gene, gene_exp, epsilon);
        w = -log2(w) + 1; % Normalized
        weight = [weight; w];            
    end
end

end
