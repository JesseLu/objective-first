function [] = int_newt_simple(fun, x, l, u, A, b, ...
                                mu_0, sigma, tau, err_tol)
% Implementation with reduced system Hessian, simple upper and lower bounds.

% Choose initial values of variables
s0 = x - l; % Slack variable for lower bound.
z0 = zeros(length(l), 1); % Dual variable for lower bound.
s1 = u - x; % Slack variable for upper bound.
z1 = zeros(length(u), 1); % Dual variable for upper bound.
y = zeros(length(b), 1); % Dual variable for equality constraint.

% Helper functions.
empty = @(x1, x2) sparse(length(x1), length(x2)); % Matrix filled with 0's.
my_diag = @(x) spdiags(x, 0, length(x), length(x)); % Sparse diagonal matrix.

% Error function.
kkt_res = @(x, s0, s1, y, z0, z1, mu) cat(1, ...
                        fun.g(x) - A' * y ...
                        - z0 + my_diag(z0./s0) * (x - l) - mu * s0.^-1 ...
                        + z1 - my_diag(z1./s1) * (u - x) + mu * s1.^-1, ...
                        A * x - b);
err = @(x, s0, s1, y, z0, z1, mu) norm(kkt_res(x, s0, s1, y, z0, z1, mu));

% System Hessian matrix.
H_sys = @(x, s0, s1, y, z0, z1, mu) ...
    [(fun.h(x) + my_diag(z0./s0) + my_diag(z1./s1)), A';
    A, empty(y, y)];

% Helper to calculate different components of p
calc_p_xy = @(p) struct(    'x', p(1:length(x)), ...
                            'y', -p(length(x) + [1:length(y)]));
calc_p_s0 = @(p, x, s0) p.x + (x - l) - s0;
calc_p_s1 = @(p, x, s1) -p.x + (u - x) - s1;
calc_p_z0 = @(p, x, s0, z0, mu) ...
            -my_diag(z0./s0) * (p.x + (x - l)) + mu * s0.^-1;
calc_p_z1 = @(p, x, s1, z1, mu) ...
            my_diag(z1./s1) * (p.x - (u - x)) + mu * s1.^-1;

% Fraction-to-boundary rule (for inequality constraint).
my_pos = @(z) (z > 0) .* z + (z <= 0) * 1; % If negative, set to 1.
f2b_rule = @(pz, z) min([1; my_pos(-tau*z./pz)]);

hist.err(1) = err(x, s0, s1, y, z0, z1, mu_0);
tic
for mu = mu_0 * sigma.^[0:100]
    for k = 1 : 100

        % Obtain search direction (p).
        p = -H_sys(x, s0, s1, y, z0, z1, mu) \ ...
            kkt_res(x, s0, s1, y, z0, z1, mu);
        p = calc_p_xy(p);
        p.s0 = calc_p_s0(p, x, s0);
        p.s1 = calc_p_s1(p, x, s1);
        p.z0 = calc_p_z0(p, x, s0, z0, mu);
        p.z1 = calc_p_z1(p, x, s1, z1, mu);

        % Compute alpha using the fraction-to-boundary rule.
        alpha_prim = f2b_rule([p.s0; p.s1], [s0; s1]);
        alpha_dual = f2b_rule([p.z0; p.z1], [z0; z1]);

        % Update variables.
        x = x + alpha_prim * p.x;
        s0 = s0 + alpha_prim * p.s0;
        s1 = s1 + alpha_prim * p.s1;
        y = y + alpha_dual * p.y;
        z0 = z0 + alpha_dual * p.z0;
        z1 = z1 + alpha_dual * p.z1;

        % Calculate error.
        hist.err(end+1) = err(x, s0, s1, y, z0, z1, mu);

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
time0 = toc
semilogy(hist.err, '.-');
xlabel('Error in KKT equations');
ylabel('iterations');
title('Interior Primal-Dual Full Newton Step Convergence');
drawnow


    %
    % Compare against cvx.
    %

tic;
path(path, genpath('../cvx'));
cvx_quiet(true);
cvx_begin
    variable x_star(length(x))
    % minimize norm(fun.A*x_star - fun.b)
    minimize fun.f_cvx(x_star)
    subject to
        A * x_star - b == 0
        x_star >= l
        x_star <= u
cvx_end
time1 = toc

fprintf('Norm difference: %1.7f%% (cvx: %e, interior_newton: %e)\n', ...
    100 * norm(x_star - x)/norm(x_star), fun.f(x_star), fun.f(x));


