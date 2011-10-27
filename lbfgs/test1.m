function [] = test1(n)
% Test lbfgs against a simple unconstrained quadratic problem.

A = spdiags(randn(n, 7), -3:3, n, n);
b = randn(n, 1);
x_star = A \ b;

% Solve using lbfgs.
fun = @(x) my_quad(A, b, x);

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


function [fval, grad] = my_quad(A, b, x)
fval = 0.5 * norm(A * x - b)^2;
grad = A' * (A * x - b);
