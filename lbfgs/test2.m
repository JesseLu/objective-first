function [] = test2(n)
% Test lbfgs against a constrained quadratic problem.

A = spdiags(randn(n, 7), -3:3, n, n);  % Use sparse matrix to speed things up.
b = randn(n, 1);
path(path, genpath('../cvx'));
cvx_quiet(true);
cvx_begin
    variable x_star(n)
    minimize norm(A * x_star - b, 2) % 2-norm
    subject to
        x_star >= 0
cvx_end

% Solve using lbfgs.
fun = @(x) my_quad_con(A, b, x, 10);

tic
res = lbfgs(fun, zeros(n, 1), ...
            'Display', 'final', ...
            'MaxIters', 1e4, ...
            'MaxFuncEvals', 1e4, ...
            'StopTol', 1e-3, ...
            'RelFuncTol', 0, ...
            'LineSearch_method', 'more-thuente');
toc


if (res.ExitFlag == 0)
    fprintf('Successfully stopped on gradient-norm termination condition.\n');
else
    error('Unsuccessful termination (exit flag %d)\n', res.ExitFlag);
end
fprintf('Error: %e (l-bfgs), %e (direct).\n', fun(res.X), fun(x_star));


function [fval, grad] = my_quad_con(A, b, x, kappa)
if any(x <= 0)
    fval = Inf;
else
    fval = 0.5 * norm(A * x - b)^2 - kappa * log(x);
end
grad = A' * (A * x - b) - kappa * x.^-1;
