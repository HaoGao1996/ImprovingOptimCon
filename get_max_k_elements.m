function [value, idx] = get_max_k_elements(A, k)
% FUCTION: get the max k elements of a matrix
% INPUT:
% A is the matrix;
% k is the top k elements;
% OUTPUT:
% value is the top k elements;
% idx is the location of top k elements;

[B, I] = sort(A(:), 1, 'descend');
value = B(1:k);
[ind_row, ind_col] = ind2sub(size(A), I(1:k));
idx = [ind_row, ind_col];

end