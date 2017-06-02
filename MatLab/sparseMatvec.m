function [b] = sparseMatvec(A, x)
%mutliplies a sparse matrix by a vector

[rows, cols, vals] = find(A);
%as close to a single-precision sparse matrix as Matlab can do
vals = single(vals);

b = zeros(size(x));

i=1;
for i = 1:nnz(A)
    b(rows(i)) = b(rows(i)) + double(vals(i))*double(x(cols(i)));
end