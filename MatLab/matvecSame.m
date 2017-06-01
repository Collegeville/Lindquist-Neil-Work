function [b] = matvecSame(A, x)
%calculates the matrix-vector multiplication
%The precision of the inputs will be used for the calculation and output

[m, n] = size(A);
b = zeros(m, 1);

for row = 1:m
    temp = 0;
    for col = 1:n
        temp = temp + A(row, col)*x(col);
    end
    b(row) = temp;
end