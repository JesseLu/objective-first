function [] = test2(n)
% Test lbfgs against a constrained quadratic problem.

    %
    % Formulate problem.
    %

A = spdiags(randn(n, 7), -3:3, n, n);  % Use sparse matrix to speed things up.
b = randn(n, 1);
l = zeros(n, 1);
u = ones(n, 1);

    
    %
    % Obtain solution using direct method.
    %

path(path, genpath('../cvx'));
cvx_quiet(true);
cvx_begin
    variable x_star(n)
    minimize norm(A * x_star - b, 2) % 2-norm
    subject to
        x_star >= l
        x_star <= u
cvx_end


    % 
    % Obtain solution using lbfgs.
    %

% Setup to run the lbfgs problem.
fun = @(x, t) my_quad_con(A, b, x, l, u, t);
t = 1;
x = mean([l, u], 2);

% Solve the homotopy problem for various barrier heights.
tic
for t = 10.^[1 : -1 : -2]
    res = lbfgs_unc(@(x) fun(x, t), x, ...
                'Display', 'off', ...
                'MaxIters', 1e4, ...
                'MaxFuncEvals', 1e4, ...
                'StopTol', 1e-3, ...
                'RelFuncTol', 0, ...
                'LineSearch_method', 'more-thuente');
    x = res.X;
    fprintf('t = %1.1e, iters = %d, time = %1.1fs\n', t, res.Iters, toc);

    % Check to see if we exited successfully or not.
    if (res.ExitFlag ~= 0)
        error('Unsuccessful termination (exit flag %d)\n', res.ExitFlag);
    end
end

% Evaluate the residuals.
fprintf('Error: %e (l-bfgs), %e (direct).\n', fun(x, 0), fun(x_star, 0));


function [fval, grad] = my_quad_con(A, b, x, l, u, t)
if any(x <= 0)
    fval = Inf;
else
    fval = 0.5 * norm(A * x - b)^2 - t * sum(log(x - l) + log(u - x));
end
grad = A' * (A * x - b) - t * ((x - l).^-1 - (u - x).^-1);
