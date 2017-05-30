function [x] = iterSolve(A, b, x0, iterationLimit)

%iterSolveSingle  Solve a linear system
%  x = iterSolve(A, b) solves for x in Ax = b
%  x = iterSolve(A, b, x0) as above, with an initial guess x0
%  x = iterSolve(A, b, x0, limit) as above, but stops after so many
%                                       iterations

if nargin < 4
    iterationLimit = 1000;
    if nargin < 3
       x0 = zeros(size(b));
    end
end

x = x0;

r = b-A*x;

c = zeros(size(b));
iteration = 0;
while (norm(r) ~= 0) && (iteration < iterationLimit)
    iteration = iteration + 1;
    c = iterSolveSingle(A, r, c, iterationLimit/20);
    x = x + double(c);
    r = b - A*x;
end
