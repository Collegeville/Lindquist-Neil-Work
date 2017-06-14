function [x] = spcg(A,b,tol,maxit,M1,M2,x0)
%SPCG   Sparse Preconditioned Conjugate Gradients Method.
%   X = SPCG(A,B) attempts to solve the system of linear equations A*X=B for
%   X. The N-by-N coefficient matrix A must be symmetric and positive
%   definite and the right hand side column vector B must have length N.
%
%   X = SPCG(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then SPCG uses the default, 1e-6.
%
%   X = SPCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then SPCG uses the default, min(N,20).
%
%   X = SPCG(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then SPCG uses the default, min(N,20).
%
%   X = SPCG(A,B,TOL,MAXIT,M) and X = SPCG(A,B,TOL,MAXIT,M1,M2) use symmetric
%   positive definite preconditioner M or M=M1*M2 and effectively solve the
%   system inv(M)*A*X = inv(M)*B for X. If M is [] then a preconditioner
%   is not applied.
%
%   X = SPCG(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess. If X0 is
%   [] then SPCG uses the default, an all zero vector.
%
%   Example:
%      n1 = 21; A = gallery('moler',n1);  b1 = A*ones(n1,1);
%      tol = 1e-6;  maxit = 15;  M = diag([10:-1:1 1 1:10]);
%      x1 = SPCG(A,b1,tol,maxit,M);
%   Or use this parameterized matrix-vector product function:
%      A = @(x,n)gallery('moler',n)*x;
%      n2 = 21; b2 = A(ones(n2,1),n2);
%      x2 = SPCG(@(x)afun(x,n2),b2,tol,maxit,M);
%
%   Class support for inputs A, M1, M2:
%      sparce matrix
%
%   Class support for inputs B,X0:
%      float: double
%
%   See also PCG

%   Copyright 1984-2013 The MathWorks, Inc.
%   Modified by Neil Lindquist, 2017

if (nargin < 2)
    error(message('MATLAB:pcg:NotEnoughInputs'));
end

%ensure b is single
b = single(b);


% Check matrix and right hand side vector inputs have appropriate sizes
[m,n] = size(A);
m = uint32(m);
n = uint32(n);
if (m ~= n)
    error(message('MATLAB:pcg:NonSquareMatrix'));
end
if ~isequal(uint32(size(b)),[m,1])
    error(message('MATLAB:pcg:RSHsizeMatchCoeffMatrix', m));
end


% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
else
    tol = double(tol)
end
warned = false;
if tol <= eps
    warning(message('MATLAB:pcg:tooSmallTolerance'));
    warned = true;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:pcg:tooBigTolerance'));
    warned = true;
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit)
    maxit = uint32(min(n,20));
else
    maxit = uint32(maxit);
end

% Check for all zero right hand side vector => all zero solution
n2b = norm(double(b));                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(n,1);                % then  solution is all zeros
    return
end


if ((nargin >= 5) && ~isempty(M1))
    existM1 = true;
    if ~isequal(uint32(size(M1)),[m,m])
        error(message('MATLAB:pcg:WrongPrecondSize', m));
    end
    if isa(M1, 'sparseSingle')
        M1 = sparse(M1);
    end
else
    existM1 = false;
end

if ((nargin >= 6) && ~isempty(M2))
    existM2 = true;
    if ~isequal(uint32(size(M2)),[m,m])
        error(message('MATLAB:pcg:WrongPrecondSize', m));
    end
    if isa(M2, 'sparseSingle')
        M2 = sparse(M2);
    end
else
    existM2 = false;
end


if ((nargin >= 7) && ~isempty(x0))
    if ~isequal(uint32(size(x0)),[n,1])
        error(message('MATLAB:pcg:WrongInitGuessSize', n));
    else
        x = single(x0);
    end
else
    x = single(zeros(n,1));
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
tolb = tol * n2b;                  % Relative tolerance
r = double(b) - sparse(A)*double(x);
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    return
end

normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = uint32(0);                          % stagnation of the method
moresteps = uint32(0);
maxmsteps = min([floor(n/50),uint32(5),n-maxit]);
maxstagsteps = uint32(3);

% loop over maxit iterations (unless convergence or failure)

for ii = uint32(1) : uint32(maxit)
    % no preconditioner
    if existM1
       y = M1 \ double(r);
       if ~all(isfinite(y))
           flag = 2;
           break;
       end
    else
        y = r;
    end
    
    if existM2
        z = single(M2 \ y);
        if ~all(isfinite(z))
            flag = 2;
            break
        end
    else % no preconditioner
        z = single(y);
    end
    
    
    rho1 = rho;
    rho = r' * double(z);
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = z;
    else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = single(double(z) + beta * double(p));
    end
    q = A*p;
    pq = double(p') * double(q);
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Check for stagnation of the method    
    if (norm(double(p))*abs(alpha) < eps*norm(double(x)))
        stag = stag + 1;
    else
        stag = uint32(0);
    end
    
    x = single(double(x) + alpha * double(p));             % form new iterate
    r = r - alpha * double(q);
    normr = norm(r);
    normr_act = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        r = double(b) - sparse(A)*double(x);
        normr_act = norm(r);
        if (normr_act <= tolb)
            flag = 0;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = uint32(0);
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:pcg:tooSmallTolerance'));
                end
                flag = 3;
                break;
            end
        end
    end
    if normr_act < normrmin      % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
    end
    if stag >= maxstagsteps
        flag = 3;
        break;
    end
end

% returned solution is first with minimal residual
if (flag ~= 0)
    r_comp = double(b) - sparse(A)*double(xmin);
    if norm(r_comp) <= normr_act
        x = xmin;
    end
end
