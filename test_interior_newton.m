function [] = test_interior_newton(n, p, err_tol)
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
H = A' * A;
fun.h = @(x) H;

fun.f_cvx = @(x) norm(A * x - b);

% Equality constraint.
A_eq = randn(p, n);
b_eq = randn(p, 1);

% Inequality constraint.
A_in = speye(n);
l = zeros(n, 1);
u = ones(n, 1);

% int_newt_red(fun, ones(n, 1), A_eq, b_eq, A_in, b_in, 1, 0.01, 0.995, err_tol);
subplot 111;
ob1_interior_newton(A, b, 0.5 * ones(n, 1), l, u, err_tol);
% ob1_interior_newton(fun, 0.5 * ones(n, 1), l, u, A_eq, b_eq, ...
%     1e-3, 0.995, 0.1, 0.5, err_tol);
