function [b] = matvecMixed(A, x)
%calculates the matrix-vector multiplication
%The output and arithmetic will be single precision

[m, n] = size(A);
b = zeros(m, 1);

for row = 1:m
    temp = single(0);
    for col = 1:n
        temp = temp + single(A(row, col))*single(x(col));
    end
    b(row) = temp;
end