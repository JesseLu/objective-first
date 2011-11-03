function [x, num_iters] = lbfgs(fun, x, tol)
n = length(x);
n_max = 5;
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
