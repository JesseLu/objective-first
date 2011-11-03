function [] = test_arrow(n, p, q)

randn('state', 1);

    %
    % Formulate problem.
    %

% Arrow matrix.
d = randn(n, 1);
A = randn(p, n);

% Low-rank matrix.
U = randn(n+p, q);
V = randn(n+p, q);

% RHS vectors.
b = randn(n+p, 1);

x = arrow_solve(d, A, U, V, b);

A_hat = fprintf('error: %e\n', ...
            norm([spdiags(d, 0, n, n), A'; A, sparse(p, p)] * x + U * (V' * x)- b));


% % Direct method.
% x_star = A_hat \ b;
