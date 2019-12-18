function [o, d, u] = get_score_of_CR(CR, DScore_gene, OCRs, cutoff_flag)
% FUNCTIONS: get the optimal influence values of a set of genes
% INPUT:
% CR is the control regions: CR{i, 1} i iteration; 
% DScore_gene is the differential score of genes;
% OCRs is the optimal control regions;
% cutoff_flag: if the elments of cutoff_flag equal to 1, these
% corresponding CRs need to participate in the calculation; Otherwise
% , vice verse.
% OUTPUT: 
% o is the matrix of optimal influence score;
% d is the matrix of desired influence score;
% u is the matrix of undesired influence score;

num_SCC = size(cutoff_flag, 1); % the number of SCCs
num_gene = size(cutoff_flag, 2); % the number of genes

d = zeros(num_SCC, num_gene);
u = zeros(num_SCC, num_gene);
all_score = sum(DScore_gene);
all_undesired = length(find(DScore_gene==0));


for i = 1:num_SCC
    for j = 1:num_gene
        if cutoff_flag(i, j)
            CR_i_j = CR{i, 1}{j, 1};
            CR_OCRs_new = unique([CR_i_j; OCRs]);
            % calculate the desired score of iteration i, gene j
            d(i, j) = sum(DScore_gene(CR_OCRs_new)) / all_score;
            % calculate the undesired score of iteration i, gene j
            u(i, j) = length(find(DScore_gene(CR_OCRs_new)==0)) / all_undesired;
        end
    end
end

o = d - u;

end