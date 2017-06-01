function [x] = iterSolveSingle(A, b, x0, limit)

%iterSolveSingle  Solve a linear system with low precision
%  x = iterSolveSingle(A, b) solves for x in Ax = b
%  x = iterSolveSingle(A, b, x0) as above, with an initial guess x0
%  x = iterSolveSingle(A, b, x0, limit) as above, but stops after so many
%                                       iterations
%
%This implementation is a Gauss-Seidel method

if nargin < 4
    limit = 100;
    if nargin < 3
       x0 = zeros(size(b));
    end
end

n = size(b);
n = n(1);

x = single(x0);
oldX = x;
A = single(A);
b = single(b);

for iteration = 1:limit    
    for i = 1:n
       x(i) = b(i);
       for j = 1:n
           if i ~= j
               x(i) = x(i) - A(i, j)*x(j);
           end
       end
       
       x(i) = x(i)/A(i, i);
    end
    
    if norm(x - oldX) < .0000001
       %seems to have stabilized, return
       break
    end
    oldX = x;
end