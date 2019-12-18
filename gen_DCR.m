function DCR = gen_DCR(A)
% FUNCTION: identify the direct control regions of genes.
% INPUT:
% A is the adjacency matrix;
% Output: 
% DCR: direct control region using maxmatching

num = length(A);
A = full(A);
DCR = cell(num, 1);

% generate a random network
A0 = zeros(num);
newlable = randperm(num);
newlable = newlable';
for p = 1:num
    for q = 1:num
        if A(p, q)
            A0(newlable(p), newlable(q)) = 1;
        end
    end
end
A0 = sparse(A0);
m0 = structural_control_max_matching(A0); 
m = zeros(num, 1);
for j = 1:num
    % revise label
    if m0(j) 
        m(find(newlable==j)) = find(newlable==m0(j));
    end
end

%% obtain new SCC and update direct control region
stems = {};
stermnodes = zeros(num, 1); % stermnodes are 1
buds = {};
budnodes = zeros(num, 1); % budnodes are 1
drivernodes = find(m==0);

% stems are disjoint
for j = 1:length(drivernodes)
    stem_j_node = drivernodes(j);
    stem_j = [];
    while(stem_j_node)
        stem_j = [stem_j; stem_j_node]; 
        stem_j_node = find(m==stem_j_node);
    end
    stermnodes(stem_j, 1) = 1;
    stems = [stems; stem_j];
end
stermnodes_list = find(stermnodes);

% budnodes and stermnodes are complementary
budnodes = 1 - stermnodes;
budnodes_list = find(budnodes);
bud_j = [];
bud_j_node = [];
bud_j_node_first = [];
while(~isempty(budnodes_list))
    if isempty(bud_j_node)
        bud_j_node = budnodes_list(1, 1);
        % delete the first node
        budnodes_list(1) = []; 

        bud_j = bud_j_node; 
        bud_j_node_first = bud_j_node;
    else     
        bud_j_node = m(bud_j_node);
        budnodes_list(find(budnodes_list==bud_j_node)) = [];

        bud_j = [bud_j; bud_j_node];

        if bud_j_node == bud_j_node_first
            buds = [buds; bud_j];
            bud_j = [];
            bud_j_node = [];
            bud_j_node_first = [];
        end
    end                      
end   
% reobtain the node list of bud
budnodes_list = find(budnodes);

% additional links (stemnodes------budnodes)
AL = [];
for p = 1:length(stermnodes_list)
    for q = 1:length(budnodes_list)
        if A(budnodes_list(q), stermnodes_list(p)) == 1
            AL = [AL; [stermnodes_list(p), budnodes_list(q)]];
        end
    end
end

% construct structural control configuration (SCC)
SCC_edge_list = [];
for j = 1:num
    if m(j)
        SCC_edge_list = [SCC_edge_list; [m(j), j]];
    end
end
SCC_edge_list = [SCC_edge_list; AL];
SCC = zeros(num);
for j = 1:length(SCC_edge_list)
    SCC(SCC_edge_list(j, 1), SCC_edge_list(j, 2)) = 1; % digraph definition
    % this is different from the definition of A
end

G = digraph(SCC);
D = distances(G);

for j = 1:num
     DCR_j = find(D(j, :)<inf);
     DCR{j, 1} = DCR_j';
end
             
end
