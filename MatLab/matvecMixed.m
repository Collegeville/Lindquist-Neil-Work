function [b] = matvecMixed(A, x)
%calculates the matrix-vector multiplication
%inputs are expected to be single precision and will be outputed as single
%precision, but the arithmetic will be done in double precision

[m, n] = size(A);
b = single(zeros(m, 1));

for row = 1:m
    temp = double(0);
    for col = 1:n
        temp = temp + double(A(row, col))*double(x(col));
    end
    b(row) = single(temp);
end