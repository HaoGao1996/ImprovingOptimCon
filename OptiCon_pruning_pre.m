function [OCNs, OCRs, OCRs_o, OCRs_idx, OCN_rate, OCNs_cutoff] = OptiCon_pruning_pre(CR, DScore_gene, DrugTarget, alpha, beta, k)
% FUNCTION: calculated the optimal control nodes (OCNs) and optimal
% control regions (OCRs) based on pruning algorithm with drug target
% preference
% INPUT:
% CR is the control regions: CR{i, 1} i iteration; 
% DScore_gene is the differential score of genes;
% DrugTarget is the set of drug target genes;
% alpha is the thereshold of grow rate of optimal influence value; 
% beta is the preferent parameter to regulate the degree of preference;
% k is the number of the initial set of start points;
% OUTPUT:
% OCNs is the set of optimal control nodes;
% OCRs is the optimal control regions;
% OCRs_o is the optimal influence score;
% OCRs_idx is the index of OCNs;
% OCN_rate is the rate of each iterations;
% OCNs_cutoff is the rate of CRs that need to be calculated.


num_SCC = length(CR); % the number of SCCs
num_gene = length(DScore_gene); % the number of genes

% construct preference matrix
DrugTarget_vec = ones(num_gene, 1);
DrugTarget_vec(DrugTarget) = beta;
DrugTarget_pre = diag(DrugTarget_vec);

OCRs_o = cell(k, 1);
OCNs = cell(k, 1);
OCRs = cell(k, 1);
OCRs_idx = cell(k, 1);
OCN_rate = cell(k, 1);
OCNs_cutoff = cell(k, 1);

% get the intial optimal infulence score
%----------------------------------------
fprintf('Calculate the initial score: \n')
%----------------------------------------

[init_o, init_d, ~] = get_score_of_CR(CR, DScore_gene, [], ones(num_SCC, num_gene));
% multiply the preference value of drug target genes
% if a gene is drug target,  its CR are more possible to 
% participate in the calculation 
init_d = init_d * DrugTarget_pre;
[top_k, idx_k] = get_max_k_elements(init_o * DrugTarget_pre, k);

for i = 1:k
    fprintf('initial ponit: %d. \n', i)
    
    OCRs_o_i = top_k(i,1); % ith  optimal influence
    OCNs_i = idx_k(i,2); % ith optimal control nodes
    OCRs_i =  CR{idx_k(i,1)}{idx_k(i,2)}; % ith optimal control regions
    OCRs_idx_i = idx_k(i, :); % ith optimal control nodes index
    OCN_rate_i = 1; % the growth rate of optimal influence score
   
    cutoff_flag = ones(num_SCC, num_gene); 
    cutoff_flag(:, OCNs_i(end)) = 0; % if a gene belongs to OCNs
    
    itr = 1;
    fprintf('iteration: %d. \n', itr)
    % Initialize the flag of cutoff
    cutoff_flag = update_cutoff_flag(cutoff_flag, init_d, OCRs_o_i(end), alpha);
    % the num of regions need to be calculated
    OCNs_cutoff_i = sum(sum(cutoff_flag)) / num_SCC / num_gene;
    
    while 1
        itr = itr+1;
        fprintf('iteration: %d. \n', itr)
        
        o_copy = OCRs_o_i(end);
        [o, ~, ~] = get_score_of_CR(CR, DScore_gene, OCRs_i, cutoff_flag);
        [OCRs_o_new, idx] = get_max_k_elements(o * DrugTarget_pre, 1);
        
        opt_rate = (OCRs_o_new - o_copy)/ o_copy;
        if opt_rate < alpha
            break;
        end
        
        OCRs_o_i = [OCRs_o_i; OCRs_o_new];
        OCNs_i = [OCNs_i; idx(1,2)];
        OCRs_i = unique([OCRs_i; CR{idx(1,1)}{idx(1,2)}]);
        OCRs_idx_i = [OCRs_idx_i; idx];
        OCN_rate_i = [OCN_rate_i; opt_rate];
        
        cutoff_flag(:, OCNs_i(end)) = 0; % if a gene belongs to OCNs
        cutoff_flag = update_cutoff_flag(cutoff_flag, init_d, OCRs_o_i(end), alpha);
        % the num of regions need to be calculated
        OCNs_cutoff_i = [OCNs_cutoff_i; sum(sum(cutoff_flag)) / num_SCC / num_gene];
    end
    
    OCRs_o{i, 1} = OCRs_o_i;
    OCNs{i, 1} = OCNs_i;
    OCRs{i, 1} = OCRs_i;
    OCRs_idx{i, 1} = OCRs_idx_i;  
    OCN_rate{i, 1} = OCN_rate_i;
    OCNs_cutoff{i, 1} = OCNs_cutoff_i;
end

end