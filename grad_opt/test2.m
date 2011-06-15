% TEST2 
%  
% Test of c-go.
% Optimization of quadratic objective with linear equality constraint.
help test2


    %
    % Specify the objective.
    %

n = 100; % Size of the problem.
A = randn(n);
b = randn(n, 1);

f = @(x) 0.5 * norm(A * x - b)^2; % Optimization objective.  
g = @(x) A' * (A * x - b); % Gradient.


    %
    % Specify the constraints.
    %

p = 5; % Number of linear constraints.

C = randn(n, p);
P = @(x) x - C * (inv(C'*C) * (C' * x)); % Projector onto space orthogonal to C.

c = @(x, grad, s) x - s * P(grad); % Constraint. This problem is unconstrained.
x0 = zeros(n, 1); % Initial starting point.


    %
    % Optimize using c-go.
    %

[x_opt, fval, ss] = opt(f, g, c, x0, 1e4);
cgo_visualize(fval, ss);


    %
    % Check against actual answer.
    %

x_check = [A'*A, C; C', zeros(p)] \ [A'*b; zeros(p,1)];
fprintf('Error in funval: %e, and var: %e\n', ...
    f(x_opt) - f(x_check(1:n)), norm(x_opt - x_check(1:n))/norm(x_check(1:n)));


