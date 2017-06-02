function [b] = sparseMatvec(A, x)
%mutliplies a sparse matrix by a vector

[m, n] = size(A);
b = zeros(n, 1);

[rows, cols, vals] = find(A);

b = zeros(n, 1);

i=1;
for i = 1:nnz(A)
    
    b(rows(i)) = b(rows(i)) + vals(i)*x(cols(i));
    
end