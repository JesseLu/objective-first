function [] = test3(n, p, err_tol)
% Test interior_newton.

    %
    % Formulate problem.
    %

% Function to minimize (quadratic, convex).
A = spdiags(randn(n, 7), -3:3, n, n);  % Use sparse matrix to speed things up.
b = randn(n, 1);
fun.f = @(x) 0.5 * norm(A * x - b)^2;
fun.g = @(x) A' * (A * x - b);
fun.h = @(x) A' * A;

% Equality constraint.
A_eq = randn(p, n);
b_eq = randn(p, 1);

% Inequality constraint.
A_in = speye(n);
b_in = zeros(n, 1);

interior_newton(fun, ones(n, 1), A_eq, b_eq, A_in, b_in, 1, 0.1, 0.995, err_tol);
