The sample code (OptiCon_pruning_pre) were partially based on the published paper,
"Hu Y, Chen C, Ding Y, et al. Optimal control nodes in disease-perturbed 
networks as targets for combination therapy[J]. Nature communications, 2019, 10(1): 2180."

In the "main_demo", we can use different type of algorithms to compare by replacing the 71 line

Case1: OptiCon (Remade by myself according to the reference paper)
[OCNs, OCRs, OCRs_o, OCRs_idx, OCN_rate] = OptiCon(CR, DScore_gene, alpha, k);

Case2: OptiCon_pruning (Speed up the original algorithm by adding pruning strategy)
[OCNs, OCRs, OCRs_o, OCRs_idx, OCN_rate, OCNs_cutoff] = OptiCon_pruning(CR, DScore_gene, alpha, k)

Case3: OptiCon_pruning_pre(Based on case2, improve the results by fusing drug-target information)
[OCNs, OCRs, OCRs_o, OCRs_idx, OCN_rate, OCNs_cutoff] = OptiCon_pruning_pre(CR, DScore_gene, DrugTarget_ID, alpha, beta, k);

This is totally free project because of my interests. 

Copyright, made by Hao Gao, 15/12/2019 
