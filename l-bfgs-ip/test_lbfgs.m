function [] = test_lbfgs(n, tol, max_iters)
% Test lbfgs against a simple unconstrained quadratic problem.

randn('state', 1);
A = spdiags(randn(n, 7), -3:3, n, n);
b = randn(n, 1);
x_star = A \ b;

% Solve using lbfgs.
fun = @(x) my_quad(A, b, x);

tic
% Custom, simple lbfgs algorithm, using lbfgs_update.
fprintf('\nCustom lbfgs algorithm\n');
[x, num_iters] = lbfgs(fun, zeros(n, 1), tol);
[f, g] = fun(x);
fprintf('%d iterations, fval = %e, ||g||/n = %e\n', num_iters, f, norm(g)/n);


function [fval, grad] = my_quad(A, b, x)
fval = 0.5 * norm(A * x - b)^2;
grad = A' * (A * x - b);
