function m = structural_control_max_matching(A)
% FUNCTION: return the matching of associated bipartite networks
% Input: A is the directed network, where aij represents that node j
% points to node i; (!!!!!!!!!!!!!!!!Note that A must be a sparse matrix)
% Output: m is the maximum matching of left nodes. If mi = 0, node i is 
% unmathed.

% Construct bipartite graph
num = length(A);
B = zeros(2 * num);
A = full(A);

for i = 1:num
    for j = 1:num
        if A(i, j)
            B(i, j+num) = 1;
            B(j+num, i) = 1;
        end
    end
end

B = sparse(B);
m0 = edmonds_maximum_cardinality_matching(B);
m = m0(1:num);
m(find(m==0)) = num;
m = m - num;

end