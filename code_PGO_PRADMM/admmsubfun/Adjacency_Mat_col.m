
function [X1,X2]=Adjacency_Mat_col(A,i)
[X1, ~, X2] = find(A(:, i));
end