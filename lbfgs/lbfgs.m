function out = lbfgs(fun, x0, varargin)

% Create parser
params = inputParser;

% Set Poblano parameters
params = poblano_params(params);

% Set parameters for this method
params.addParamValue('M',5,@(x) x > 0);

% Parse input
params.parse(varargin{:});

%% Check input arguments
if (nargin == 1) && isequal(fun,'defaults') && (nargout == 1)
    out = params.Results;
    return;
elseif (nargin < 2)
    error('Error: invalid input arguments');
end

%% Initialize

xk = x0;
[fk,gk] = feval(fun,xk);
out = poblano_out(xk,fk,gk,1,params);

%% Main loop
while out.ExitFlag == -1

    if out.Iters == 0
        % Initialize quantities before first iteration
        pk = -gk;
        ak = 1.0;
        S = [];
        Y = [];
    else
        % Precompute quantites used in this iteration
        sk = xk - xkold;
        yk = gk - gkold;
        skyk = sk'*yk;
        ykyk = yk'*yk;
        rhok = 1 / skyk;
        gamma = skyk/ykyk;

        % Check curvature condition.
        if (skyk <= 0)
            error('Curvature condition broken (sk^T yk = %e)!', skyk);
        end

        % Use information from last M iterations only
        if out.Iters <= params.Results.M
            S = [S sk];
            Y = [Y yk];
        else
            S = [S(:,2:end) sk];
            Y = [Y(:,2:end) yk];
        end
        % Adjust M to available number of iterations 
        m = size(S,2);

        % Outer-product (or compact) representation.
        n = size(S, 1);
        gamma = gamma.^-1;
        D = diag(dot(S, Y, 1));
        L = S' * Y;
        L = L - triu(L); % Only keep elements below the diagonal.
        % L = L - diag(diag(L));
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
    [xk,fk,gk,ak,lsinfo,nfev] = poblano_linesearch(fun,xk,fk,gk,ak,pk,params.Results);
    f = @(t) fun(xkold + t * pk);
    t = backtrack_linesearch(f, 1, f(0), gkold'*pk, 0.1, 0.5);
    xk = xkold + t * pk;
    [fk, gk] = fun(xk);
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

