function [] = test1(n, tol, max_iters)
% Test lbfgs against a simple unconstrained quadratic problem.

randn('state', 1);
A = spdiags(randn(n, 7), -3:3, n, n);
b = randn(n, 1);
x_star = A \ b;

% Solve using lbfgs.
fun = @(x) my_quad(A, b, x);

trials = {  {'lbfgs_compact', 'Compact with MT linesearch'}, ...
            {'lbfgs_outer', 'Uses lbfgs_update with MT linesearch'}, ...
            {'lbfgs_compact_bt', 'Compact with BT linesearch'}, ...
            {'lbfgs_outer_bt', 'Uses lbfgs_update with BT linesearch'}};

for k = 1 : length(trials)
    fprintf(['\n', trials{k}{2}, '\n']);
    tstart = tic;
    res = feval(trials{k}{1}, fun, zeros(n, 1), ...
                'Display', 'off', ...
                'MaxIters', max_iters, ...
                'MaxFuncEvals', 1e4, ...
                'StopTol', tol, ...
                'RelFuncTol', 0, ...
                'LineSearch_method', 'more-thuente');
    fprintf('iters: %d, fevals: %d, fval: %1.3e, err: %1.3e, time: %1.3f s\n', ...
        res.Iters, res.FuncEvals, res.F, norm(res.G)/n, toc(tstart));
end

tic
% Custom, simple lbfgs algorithm, using lbfgs_update.
fprintf('\nCustom lbfgs algorithm\n');
n_max = 5;
x = zeros(n, 1);
[f, g] = feval(fun, x);
h = [];
num_iters = 0;
while norm(g)/n > tol % Termination condition.
    if isempty(h) % Initialize.
        [delta, M, W, h] = lbfgs_update(x, g, n_max, h);
        p = -g; % Just go in steepest-descent direction.
    else
        [delta, M, W, h] = lbfgs_update(x, g, n_max, h);
        p = arrow_solve((delta)*ones(n,1), ones(0,n), -W*M, W, -g);
    end
    
    f_t = @(t) fun(x + t * p);
    t = backtrack_linesearch(f_t, 1, f_t(0), (g' * p), 0.1, 0.5);

    x_prev = x;
    g_prev = g;

    x = x_prev + t * p;
    [f, g] = fun(x);
    num_iters = num_iters + 1;
end
fprintf('%d iterations, fval = %e, ||g||/n = %e\n', num_iters, f, norm(g)/n);


function [fval, grad] = my_quad(A, b, x)
fval = 0.5 * norm(A * x - b)^2;
grad = A' * (A * x - b);
