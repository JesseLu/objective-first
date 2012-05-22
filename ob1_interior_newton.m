%% Interior-point Newton algorithm for simple bounds.
function [x, err] = ob1_interior_newton(fun, x, l, u, A, b, ...
                                sigma, tau, alpha, beta, err_tol)


        %
        % Set up variables and helper functions.
        %

    % Choose initial values of variables
    s0 = x - l; % Slack variable for lower bound.
    z0 = ones(length(l), 1); % Dual variable for lower bound.
    s1 = u - x; % Slack variable for upper bound.
    z1 = ones(length(u), 1); % Dual variable for upper bound.
    y = zeros(length(b), 1); % Dual variable for equality constraint.

    % Helper functions.
    my_diag = @(x) spdiags(x, 0, length(x), length(x)); % Sparse diagonal matrix.

    % RHS of KKT equation.
    kkt_res = @(x, s0, s1, y, z0, z1, mu) cat(1, ...
                            fun.g(x) - A' * y ...
                            - z0 + (z0./s0) .* (x - l) - mu * s0.^-1 ...
                            + z1 - (z1./s1) .* (u - x) + mu * s1.^-1, ...
                            A * x - b);

    % System Hessian matrix.
    empty = sparse(length(y), length(y)); % Matrix filled with 0's.
    H_sys = @(x, s0, s1, y, z0, z1, mu) ...
        [(fun.h(x) + my_diag(z0./s0 + z1./s1)), A';
        A, empty];

    % Merit function.
    phi = @(x, s0, s1, y, z0, z1, mu, alpha_prim, alpha_dual, p, t) ...
        norm(kkt_res(   x + t * alpha_prim * p.x, ...
                        s0 + t * alpha_prim * p.s0, ...
                        s1 + t * alpha_prim * p.s1, ...
                        y + 1 * alpha_dual * p.y, ...
                        z0 + 1 * alpha_dual * p.z0, ...
                        z1 + 1 * alpha_dual * p.z1, mu));

    % Error function.
    err = @(x, s0, s1, y, z0, z1, mu) max(cat(1, ...
                            norm(fun.g(x) - A' * y - z0 + z1), ...
                            norm(s0 .* z0 - mu), ...
                            norm(s1 .* z1 - mu), ...
                            norm(A * x - b), ...
                            norm(x - l - s0), ...
                            norm(u - x - s1)));

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


        %
        % Optimize!
        %

    hist.err(1) = err(x, s0, s1, y, z0, z1, 0);
    hist.err_mu(1) = err(x, s0, s1, y, z0, z1, sigma * hist.err(1));
    hist.t(1) = nan;
    tic

    for k = 1:100
        mu = sigma * hist.err(end); % Dynamically set mu.

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

        % Perform a backtracking (Armijo) line search.
        t = simple_backtrack(@(t) ...
            phi(x, s0, s1, y, z0, z1, mu, alpha_prim, alpha_dual, p, t), ...
            alpha_prim, alpha, beta);

        % Update variables.
        x = x + t * alpha_prim * p.x;
        s0 = s0 + t * alpha_prim * p.s0;
        s1 = s1 + t * alpha_prim * p.s1;
        y = y + 1 * alpha_dual * p.y;
        z0 = z0 + 1 * alpha_dual * p.z0;
        z1 = z1 + 1 * alpha_dual * p.z1;

        % Calculate error.
        hist.err(end+1) = err(x, s0, s1, y, z0, z1, 0);
        hist.err_mu(end+1) = err(x, s0, s1, y, z0, z1, mu);
        hist.t(end+1) = t;

        % Test termination condition.
        if (hist.err(end) <= err_tol)
            break
        end
    end

    time0 = toc;
    hist.err(end)
    % Plot results.
    semilogy(0:length(hist.err)-1, [hist.err; hist.err_mu; hist.t]', '.-');
    xlabel('Error in KKT equations');
    ylabel('iterations');
    title('Interior Primal-Dual Full Newton Step Convergence');
    legend({'global error', 'local error', 'step size'}, -1);
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
    time1 = toc;

    fprintf('Percent error: %1.7f%% (compared against cvx result)\n',  ...
        100 * norm(x_star - x)/norm(x_star));
    fprintf('Interior newton, fval: %e, time: %1.2f s\n', fun.f(x), time0);
    fprintf('cvx, fval: %e, time: %1.2f s\n', fun.f(x_star), time1);


function [t] = simple_backtrack(f, t, alpha, beta);
% Backtracking line search on one-dimensional function f.
    % t = 1;
    f0 = f(0);
    y = []; % For debugging purposes.
    x = [];
    while f(t) > (1 - alpha * t) * f(0)
        y(end+1) = f(t);
        x(end+1) = t;
        t = beta * t;
        if (t <= 1e-6) % Just try to get some improvement.
            warning('Setting alpha to 0');
            alpha = 0;
        end
        if (t <= eps) % Not a descent direction.
            semilogx(x, y - f0, '.-');
            drawnow
            error('Backtracking line-search failed.');
            break
        end
    end

