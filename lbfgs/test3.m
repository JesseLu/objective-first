function [] = test3(n, p, err_tol)
% Test interior_newton.

randn('state', 1);

    %
    % Formulate problem.
    %

% Function to minimize (quadratic, convex).
A = spdiags(randn(n, 7), -3:3, n, n);  % Use sparse matrix to speed things up.
b = randn(n, 1);
fun.f = @(x) 0.5 * norm(A * x - b)^2;
fun.g = @(x) A' * (A * x - b);
fun.h = @(x) A' * A;

fun.f_cvx = @(x) norm(A * x - b);

% Equality constraint.
A_eq = randn(p, n);
b_eq = randn(p, 1);

% Inequality constraint.
A_in = speye(n);
b_in = zeros(n, 1);

subplot 121;
interior_newton(fun, ones(n, 1), A_eq, b_eq, A_in, b_in, 1, 0.01, 0.995, err_tol);
fprintf('---\n')
subplot 122;
int_newt_red(fun, ones(n, 1), A_eq, b_eq, A_in, b_in, 1, 0.01, 0.995, err_tol);
