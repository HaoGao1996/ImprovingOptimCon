function ICR = gen_ICR(A, gene_exp, DScore_gene, lambda, epsilon)
% FUNCTION: generate the indirect control regions of genes.
% Input:
% A is the adjacency matrix;
% gene_exp is the expression data of genes in GRN;
% DScore_gene is the processed differential expression genes;
% lambda is the threshold of ICV
% epsilon is the value correpondings to a p-value of 0.05 based on the
% theoretical distribution of the PCC.
% Output: 
% ICR: indirect control region of genes

num = length(A); % num of genes
ICR = cell(num, 1); % store the results

for i = 1:num
    %---------------------------------------------
    fprintf('Genrate ICR of gene.%d \n', i)
    %---------------------------------------------
    if DScore_gene(i)
        der_gene = i;
        % construct weighted downstream subnetwork
        [weightednet, downnodes] = get_weightednet(A, gene_exp, der_gene, epsilon);
        
        ICR{i, 1} = get_ICR_der_gene(weightednet, downnodes, der_gene, lambda);
    end
end
           
end
