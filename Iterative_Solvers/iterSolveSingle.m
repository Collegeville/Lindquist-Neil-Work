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

x = single(x0);
A = single(A);
b = single(b);

invL = inv(tril(A));
T = -invL * triu(A, 1);
C = invL*b;
iteration = 0;
while iteration < limit
    iteration = iteration +1;
    newX = T*x+C;
    if norm(x - newX) < .000001
       %seems to have stabilized, return
       x = newX;
       break
    end
    x = newX;
end