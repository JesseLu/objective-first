% TEST 1
%  
%  Simple test of c-go. Gradient optimization of a convex function.


    %
    % Specify the problem.
    %

n = 100; % Size of the problem.
A = randn(n);
b = randn(n, 1);

f = @(x) 0.5 * norm(A * x - b)^2; % Optimization objective.  
g = @(x) A' * (A * x - b); % Gradient.
c = @(x, grad, s) x - s * grad; % Constraint. This problem is unconstrained.
x = zeros(n, 1); % Initial starting point.


    %
    % Optimize using c-go.
    %

x_opt = opt(f, g, c, x);


    %
    % Evaluate how close we got to the actual answer.
    %
