function cutoff_flag = update_cutoff_flag(cutoff_flag, init_d, o_k, alpha)
% FUNCTION: updata the cutoff_flag matrix through init_d and pre-optimal
% influence
% INPUT: 
% cutoff_flag: if the elments of cutoff_flag equal to 1, these
% corresponding CRs need to participate in the calculation; Otherwise
% , vice verse.
% init_d is the initial desired influence score matrix. 
% o_k is the optimal influence score of the kth iteration.
% alpha is the terminal condition of the growth rate of the optinmal 
% influence score.
% OUTPUT:
% cutoff_flag

num_SCC = size(cutoff_flag, 1); % the number of SCCs
num_gene = size(cutoff_flag, 2); % the number of genes

for i = 1:num_SCC
    for j = 1:num_gene
        % update the flag of CR if its cutoff_flag value is 1 and it
        % satisfies the inequality.
        if cutoff_flag(i, j) && init_d(i, j) / o_k <= alpha
            cutoff_flag(i, j) = 0;          
        end
    end
end

end