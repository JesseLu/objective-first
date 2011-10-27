function [] = test2(n)
% Test lbfgs against a constrained quadratic problem.

    %
    % Formulate problem.
    %

A = spdiags(randn(n, 7), -3:3, n, n);  % Use sparse matrix to speed things up.
b = randn(n, 1);
l = zeros(n, 1);
u = ones(n, 1);

    
    %
    % Obtain solution using direct method.
    %

path(path, genpath('../cvx'));
cvx_quiet(true);
cvx_begin
    variable x_star(n)
    minimize norm(A * x_star - b, 2) % 2-norm
    subject to
        x_star >= l
        x_star <= u
cvx_end


    % 
    % Obtain solution using lbfgs.
    %

% Setup to run the lbfgs problem.
fun = @(x, t) my_quad_con(A, b, x, l, u, t);
t = 1;
x = mean([l, u], 2);

% % Verify gradient.
% out = gradientcheck(@(x) fun(x, 1), x);
% out.NormGradientDiffs

% Solve the homotopy problem for various barrier heights.
hist.fevals = [];
hist.fun = [];
tic;
for t = 10.^-[1:0.5:4]
    res = lbfgs_con(@(x) fun(x, t), x, my_H0(A, b, l, u, t), ...
                'Display', 'off', ...
                'MaxIters', 1e4, ...
                'MaxFuncEvals', 1e4, ...
                'StopTol', t, ...
                'RelFuncTol', 0, ...
                'LineSearch_maxfev', 1e2, ...
                'LineSearch_method', 'more-thuente');
    x = res.X;
    hist.fevals(end+1) = res.FuncEvals;
    hist.fun(end+1) = fun(x, 0);

    fprintf('t = %1.1e, fevals = %d, time = %1.1fs\n', t, res.FuncEvals, toc);

    % Check to see if we exited successfully or not.
    if (res.ExitFlag ~= 0)
        error('Unsuccessful termination (exit flag %d)\n', res.ExitFlag);
    end
end

% Evaluate the residuals.
fprintf('Error: %e (l-bfgs), %e (direct).\n', fun(x, 0), fun(x_star, 0));
subplot 211; semilogy(cumsum(hist.fevals), hist.fun - fun(x_star, 0), '.-');
subplot 212; plot([x_star, x], '.-');


function [fval, grad] = my_quad_con(A, b, x, l, u, t)
if any(x <= l) | any(x >= u)
    fval = Inf;
else
    fval = 0.5 * norm(A * x - b)^2 - t * sum(log(x - l) + log(u - x));
end
grad = A' * (A * x - b) - t * ((x - l).^-1 - (u - x).^-1);

function [H0] = my_H0(A, b, l, u, t)
H0 = @(x, q) q ./ (t * ((x - l).^-2  + (u - x).^-2));
