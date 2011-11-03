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
[x, num_iters] = lbfgs(fun, zeros(n, 1), tol);
[f, g] = fun(x);
fprintf('%d iterations, fval = %e, ||g||/n = %e\n', num_iters, f, norm(g)/n);


function [fval, grad] = my_quad(A, b, x)
fval = 0.5 * norm(A * x - b)^2;
grad = A' * (A * x - b);
