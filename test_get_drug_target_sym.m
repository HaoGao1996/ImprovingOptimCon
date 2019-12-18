clear
clc
close

load LungData
% load final_OCN
load final_OCN_pruning_pre

OCNs_DT = OCNs{1,1};
[isDT, ~] = ismember(OCNs_DT, DrugTarget_ID);
OCNs_DT = OCNs_DT .* isDT;
OCNs_DT (OCNs_DT ==0 ) = [];
OCNs_DT_ID_en = GeneID(OCNs_DT);
[~, loc] = ismember(OCNs_DT_ID_en, symbol2entrez_Integ(:, 2));
OCNs_DT_ID_sym = symbol2entrez_Integ(loc, 1)
