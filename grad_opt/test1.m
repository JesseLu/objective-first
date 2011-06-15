% TEST1
%  
% Simple test of c-go. Gradient optimization of a convex function.
help test1

    %
    % Specify the problem.
    %

n = 100; % Size of the problem.
A = randn(n);
b = randn(n, 1);

f = @(x) 0.5 * norm(A * x - b)^2; % Optimization objective.  
g = @(x) A' * (A * x - b); % Gradient.
c = @(x, grad, s) x - s * grad; % Constraint. This problem is unconstrained.
x0 = zeros(n, 1); % Initial starting point.


    %
    % Optimize using c-go.
    %

[x_opt, fval, ss] = opt(f, g, c, x0, 5e4);
cgo_visualize(fval, ss);


    %
    % Check against actual answer.
    %

x_check = A \ b;
fprintf('Error in funval: %e, and var: %e\n', ...
    f(x_opt) - f(x_check), norm(x_opt - x_check)/norm(x_check));


