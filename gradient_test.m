function [avg_err, err] = gradient_test(f, g, x0, v_mag, num_iters)

f0 = f(x0);
g0 = g(x0);
n = length(g0);

for k = 1 : num_iters
    % Get random vector v.
    v = randn(n, 1);
    v = v_mag.^-1 * (v ./ norm(v));

    % Compute error.
    err(k) = norm((f0 - f(x0 + v)) - real(g0' * v));
end

% Average error.
avg_err = mean(err);
