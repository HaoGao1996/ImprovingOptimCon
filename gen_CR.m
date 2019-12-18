function CR = gen_CR(A, ICR, DScore_gene, num_itr, GeneID)
% Function: gent the control regions of genes
% INPUT:
% ICR is the indirected control region of genes;
% num_itr is the number of iteration to find SCC;
% GeneID is the genes in GRN;
% OUTPUT:
% CR is the set of control regions compute by different SCC.

num_genes = length(GeneID);
CR = cell(num_itr, 1);

for i = 1:num_itr
    %---------------------------------------------
    fprintf('Genrate CR of gene.%d \n', i)
    %---------------------------------------------
    CR_i = cell(num_genes, 1);
    DCR = gen_DCR(A);  
    
    for j = 1:num_genes
        CR_i_j = DCR{j, 1};
        for k = 1:length(DCR{j, 1})       
            if DScore_gene(DCR{j, 1}(k)) ~= 0
                CR_i_j = [CR_i_j; ICR{DCR{j, 1}(k), 1}];
            end
        end
        CR_i{j, 1} = unique(CR_i_j);
    end
    
    CR{i, 1} = CR_i;   
end

end