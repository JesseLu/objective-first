function [] = interior_newton(fun, x, A_eq, b_eq, A_in, b_in, ...
                                mu_0, sigma, tau, err_tol)

% Choose initial values of variables
s = x-b_in;
y = ones(length(b_eq), 1); % Dual variable for equality constraint.
z = ones(length(b_in), 1); % Dual variable for inequality constraint.

% Error function.
kkt_res = @(x, s, y, z, mu) cat(1, ...
                                fun.g(x) - A_eq' * y - A_in' * z, ...
                                z - mu * s.^-1, ...
                                A_eq * x - b_eq, ...
                                A_in * x - b_in - s);
err = @(x, s, y, z, mu) norm(kkt_res(x, s, y, z, mu));

% Helper functions (for creating system Hessian).
empty = @(x1, x2) sparse(length(x1), length(x2)); % Matrix filled with 0's.
my_diag = @(x) spdiags(x, 0, length(x), length(x)); % Sparse diagonal matrix.

% System Hessian matrix.
H_sys = @(x, s, y, z, mu) ...
    [fun.h(x), empty(x, z), A_eq', A_in'; ...
    empty(s, x), my_diag(z./s), empty(s, y), -speye(length(s)); ...
    A_eq, empty(y, [s; y; z]); ...
    A_in, -speye(length(s)), empty(z, [y; z])];

% Helper function to split up p.
split_p = @(p) struct(  'x', p(1:length(x)), ...
                        's', p(length(x) + [1:length(s)]), ...
                        'y', -p(length([x; s]) + [1:length(y)]), ...
                        'z', -p(length([x; s; y]) + [1:length(z)]));

% Fraction-to-boundary rule (for inequality constraint).
my_pos = @(z) (z > 0) .* z + (z <= 0) * 1; % If negative, set to 1.
f2b_rule = @(pz, z) min([1; my_pos(-tau*z./pz)]);

hist.err(1) = err(x, s, y, z, mu_0);

for mu = mu_0 * sigma.^[0:100]
    for k = 1 : 100

        % Obtain search direction (p).
        p = -H_sys(x, s, y, z, mu) \ kkt_res(x, s, y, z, mu);
        p = split_p(p);

        % Compute alpha using the fraction-to-boundary rule.
        alpha_prim = f2b_rule(p.s, s);
        alpha_dual = f2b_rule(p.z, z);

        % Update variables.
        x = x + alpha_prim * p.x;
        s = s + alpha_prim * p.s;
        y = y + alpha_dual * p.y;
        z = z + alpha_dual * p.z;

        % Calculate error.
        hist.err(end+1) = err(x, s, y, z, mu);

        % Test inner loop termination condition.
        if (hist.err(end) <= mu)
            break
        end
    end

    % Test outer loop termination condition.
    if (mu <= err_tol)
        break
    end
end

semilogy(hist.err, '.-');
xlabel('Error in KKT equations');
ylabel('iterations');
title('Interior Primal-Dual Full Newton Step Convergence');


    %
    % Compare against cvx.
    %

path(path, genpath('../cvx'));
cvx_quiet(true);
cvx_begin
    variable x_star(length(x))
    % minimize norm(fun.A*x_star - fun.b)
    minimize fun.f_cvx(x_star)
    subject to
        A_eq * x_star - b_eq == 0
        A_in * x_star - b_in >= 0
cvx_end

fprintf('Percent difference: %1.7f%% (cvx: %e, interior_newton: %e)\n', ...
    100 * (fun.f(x_star) - fun.f(x))/fun.f(x_star), fun.f(x_star), fun.f(x));
