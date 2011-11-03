function out = lbfgs_compact(FUN,x0,varargin)
%LBFGS   Limited-memory BFGS minimization (vector-based).
%
%  OUT = LBFGS(FUN,X0) minimizes FUN starting at the point X0 using L-BFGS.
%  FUN is a handle for a function that takes a single vector input and
%  returns two arguments --- the scalar function value and the
%  vector-valued gradient. See POBLANO_OUT for details of the output
%  parameters.
%
%  OUT = LBFGS(FUN,X0,'param',value,...) specifies a parameters and its
%  value. See POBLANO_PARAMS for further details on standard parameters.
%  Additionally, LBFGS requires
%
%   'M' - Limited memory parameter {5}.
%
%  PARAMS = LBFGS('defaults') returns a structure containing the 
%  default parameters for the particular Poblano method. 
%
%
%  Examples 
%  
%  Suppose the function and gradient of the objective function are
%  specified in an mfile named example1.m:
%
%    function [f,g]=example1(x,a)
%    if nargin < 2, a = 1; end
%    f = sin(a*x);
%    g = a*cos(a*x);
%
%  We can call the optimization method (using its default
%  parameters) using the command:
%
%    out = lbfgs(@(x) example1(x,3), pi/4);
%
%  To change a parameter, we can specify a param/value input pair
%  as follows:
%
%    out = lbfgs(@(x) example1(x,3), pi/4, 'Display', 'final');
%
%  Alternatively, we can use a structure to define the parameters:
%  
%    params.MaxIters = 2;
%    out = lbfgs(@(x) example1(x,3), pi/4, params);
%
%  See also POBLANO_OUT, POBLANO_PARAMS, FUNCTION_HANDLE.
%
%MATLAB Poblano Toolbox.
%Copyright 2009, Sandia Corporation.

%% Parse parameters

% Create parser
params = inputParser;

% Set Poblano parameters
params = poblano_params(params);

% Set parameters for this method
params.addParamValue('M',5,@(x) x > 0);

% Parse input
params.parse(varargin{:});

%% Check input arguments
if (nargin == 1) && isequal(FUN,'defaults') && (nargout == 1)
    out = params.Results;
    return;
elseif (nargin < 2)
    error('Error: invalid input arguments');
end

%% Initialize

xk = x0;
[fk,gk] = feval(FUN,xk);
out = poblano_out(xk,fk,gk,1,params);

%% Main loop
while out.ExitFlag == -1

    if out.Iters == 0
        % Initialize quantities before first iteration
        pk = -gk;
        ak = 1.0;
        S = [];
        Y = [];
        rho = [];
    else
        % Precompute quantites used in this iteration
        sk = xk - xkold;
        yk = gk - gkold;
        skyk = yk'*sk;
        ykyk = yk'*yk;
        rhok = 1 / skyk;
        gamma = skyk/ykyk;

        % Use information from last M iterations only
        if out.Iters <= params.Results.M
            S = [S sk];
            Y = [Y yk];
            rho = [rho rhok];
        else
            S = [S(:,2:end) sk];
            Y = [Y(:,2:end) yk];
            rho = [rho(2:end) rhok];
        end
        % Adjust M to available number of iterations 
        m = size(S,2);

%         % L-BFGS two-loop recursion
%         q = gk;        
%         for i = 1:m
%             alpha(i) = rho(i)*S(:,i)'*q;
%             q = q - alpha(i)*Y(:,i);
%         end
%         r = gamma*q;
%         for i = m:-1:1
%             beta = rho(i)*Y(:,i)'*r;
%             r = r + (alpha(i) - beta)*S(:,i);
%         end

        % Outer-product (or compact) representation.
        n = size(S, 1);
        gamma = gamma.^-1;
        D = diag(dot(S, Y, 1));
        L = tril(S' * Y);
        L = L - diag(diag(L));
        M = inv([gamma*S'*S, L; L', -D]);
        W = [gamma*S, Y];
        Ainv = (1/gamma) * speye(n);
        r = inv_lemma(Ainv, -W, M*W', gk);
        
        % r contains H_k * g_k (Hessian approximation at iteration k times
        % the gradient at iteration k
        pk = -r;
        
    end
    xkold = xk;
    gkold = gk;

    % Compute step length
    [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch(FUN,xk,fk,gk,ak,pk,params.Results);
    f = @(t) FUN(xkold + t * pk);
    t = backtrack_linesearch(f, 1, f(0), gkold'*pk, 0.1, 0.5);
    xk = xkold + t * pk;
    [fk, gk] = FUN(xk);
    if (lsinfo ~= 1) && strcmp(params.Results.Display, 'iter')
        fprintf(1,[mfilename,': line search warning = %d\n'],lsinfo);
    end
    
    % Update counts, check exit conditions, etc.
    out = poblano_out(xk,fk,gk,nfev,params,out);    
end

function [x] = inv_lemma(Ainv, B, C, b)
% Matrix inversion lemma (also known as Sherman-Woodbury-Morrison formula).
% An efficient way of computing (A + BC)x = b, given A^-1.
% Reference: Appendix C.4.3, Convex Optimization, Boyd and Vandenberghe
z = Ainv * b;
E = speye(size(C, 1)) + C * Ainv * B;
w = E \ (C * z);
x = z - Ainv * B * w;

