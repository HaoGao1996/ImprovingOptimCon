function w = comp_w(u, v, j, gene_exp, epsilon)
% FUNCTION: get the weighted network of downstream nodes
% Input:
% edge <u, v> is in the downstream subnetwork of j;
% j is the deregulated gene;
% gene_exp is the expression data of genes in GRN;
% epsilon is the value correpondings to a p-value of 0.05 based on the
% theoretical distribution of the PCC.
% Output:
% w is the weight of edge <u, v> of j

gj=gene_exp(j,:); 
gu=gene_exp(u,:);
gv=gene_exp(v,:);

% gj, gu and gv may be all 0;
if ~sum(gj)
    PCC1=0;
    PCC2=0;
else  
 
    if ~sum(gu)
        PCC1 = 0;
    else
        R1 = corrcoef(gj',gu');
        PCC1 = abs(R1(1,2));
    end
    
    if ~sum(gv)
        PCC2 = 0;
    else
        R2 = corrcoef(gj',gv');
        PCC2 = abs(R2(1,2));
    end
end

w = (PCC1 + PCC2) / 2;
if w < epsilon
    w = epsilon;
end

end

