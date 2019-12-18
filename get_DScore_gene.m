function [DScore_gene, count] = get_DScore_gene(GeneID, DScore_ID, DScore_all)
% Function: Get the expression matrix of genes in GRN
% Input:
% GeneID: the gene ID in GRN;
% DScore_ID: all genes ID
% DScore_all: the associated DScore of all genes
% Output:
% DScore_gene: the DScore of genes in GRN;
% count: the number of genes that have DScore;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Genes without expression data

DScore_gene = zeros(length(GeneID), 1);
count=0;    

for i=1:length(GeneID)
    idx = find(strcmp(GeneID(i),DScore_ID));
    if ~isempty(idx)  
        if DScore_all(idx) < 0.05
            DScore_gene(i, :) = -log10(DScore_all(idx));  
            count = count + 1;
        end
    end
end

end

