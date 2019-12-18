clear
close
clc

% if there is no toolkit matlab_bgl in the folder, please download it 
% from https://www.mathworks.com/matlabcentral/fileexchange/10922-matlabbgl
% addpath('E:\Code\Matlab\OptiCon\matlab_bgl\');

%% load data and process data
%---------------------------------------------
fprintf('Load data: \n')
%---------------------------------------------
% Read network from edgelist
inputdatadirection = 'E:\Code\Matlab\OptiCon_Pre\data\';
[a1, a2] = textread([inputdatadirection 'Network_original.txt'], '%s%s');
GeneName = unique([a1; a2]);
[~, loc_a1] = ismember(a1, GeneName);
[~, loc_a2] = ismember(a2, GeneName);
GeneEdgeList = [loc_a1, loc_a2];
adjMatrix = zeros(length(GeneName)); % generate adjacent matrix
for i = 1:length(a1)
    adjMatrix(loc_a2(i, 1), loc_a1(i, 1)) = 1; % a1----->a2
end
adjMatrix = sparse(adjMatrix);

% Read original gene expression data (FPKM)
fpkm = bioma.data.DataMatrix('File', [inputdatadirection, 'GeneExpression.txt']);
gene_exp = get_gene_exp(GeneID,fpkm);

% Read the adjusted p-value of differential expressed genes
[DScore_ID, DScore_all] = textread([inputdatadirection, 'DScore.txt'],'%s%f');  
% Map the genes in Gene regulatory network
DScore_gene = get_DScore_gene(GeneID, DScore_ID, DScore_all);

% Read the recurrent mutated genes
CancerGene_en = importdata([inputdatadirection, 'CancerCensus.txt']);
[~, CancerGene_ID] = ismember(CancerGene_en, GeneID);
CancerGene_ID(CancerGene_ID == 0) = [];
% Read Drug-target

DrugTarget_symbol = importdata([inputdatadirection, 'ApprovedDrugTarget.txt']);
[~, DrugTarget_en_ID] = ismember(DrugTarget_symbol, symbol2entrez_Integ(:, 1));
DrugTarget_en_ID(DrugTarget_en_ID == 0) = [];

DrugTarget_en = symbol2entrez_Integ(DrugTarget_en_ID, 2);
[~, DrugTarget_ID] = ismember(DrugTarget_en, GeneID);
DrugTarget_ID(DrugTarget_ID == 0) = [];
clearvars -except LungData adjMatrix DrugTarget_ID CancerGene_ID DScore_gene gene_exp GeneID integ_RN symbol2entrez_Integ

%% Generate the indirect control regions
%---------------------------------------------
fprintf('Genrate ICR: \n')
%---------------------------------------------
lambda = 0.3; % parameter 1
epsilon = 0.18; % parameter 2
ICR = gen_ICR(adjMatrix, gene_exp, DScore_gene, lambda, epsilon);

%% Generate Control regions according to random sampling method
num_itr = 1000; % times of iterations
%---------------------------------------------
fprintf('Genrate CR: \n')
%---------------------------------------------
CR = gen_CR(adjMatrix, ICR, DScore_gene, num_itr, GeneID);

%% Identify the OCNs with preference
% k = ceil(0.0001 * length(CR) * length(DScore_gene));
k = 1; 
alpha = 0.05;
beta = 1.07; % beta need to be no less than 1;
tic
[OCNs, OCRs, OCRs_o, OCRs_idx, OCN_rate, OCNs_cutoff] = OptiCon_pruning_pre(CR, DScore_gene, DrugTarget_ID, alpha, beta, k);
toc

save final_OCN_pruning_pre OCNs OCRs OCRs_o OCRs_idx OCN_rate OCNs_cutoff


