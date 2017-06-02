function [b] = sparseMatvec(A, x)
%mutliplies a sparse matrix by a vector

[rows, cols, vals] = find(A);

b = zeros(size(x));

for i = 1:nnz(A)
    b(rows(i)) = b(rows(i)) + vals(i)*double(x(cols(i)));
end