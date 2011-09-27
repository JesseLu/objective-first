function [epsilon] =  ob1_demo1(omega, epsilon, eps_lims, out_dir, ...
    mode_nums, num_iters, options)

dims = size(epsilon);
N = prod(dims);
path(path, genpath('~/lumos/cvx'));


    %
    % Form matrices.
    %

[A, S] = ob1_matrices(dims, dims-4);


    % 
    % Get initial value of p.
    %

% Start with maximum epsilon within active box.
p0 = 1 ./ (epsilon(:) - (S.p * S.p') * (epsilon(:) - options(1)));

    
    % 
    % Get initial value of x.
    %

eps = A{2} * epsilon(:);
eps = struct('x', reshape(eps(1:N), dims), 'y', reshape(eps(N+1:2*N), dims));

mode = {ob1_wgmode(omega, eps, 'x-', 'in', mode_nums(1)), ...
    ob1_wgmode(omega, eps, out_dir, 'out', mode_nums(2))};
x0 = mode{1} + mode{2};


    %
    % Optimize.
    %

% Helper function for sparse diagonal matrix.
D_ = @(z) sparse(1:length(z), 1:length(z), z, length(z), length(z));

% Physics residual calculation.
phys_res = @(x, p) ...
    norm(S.r' * (A{1} * D_(A{2} * p) * A{3} * x - omega^2 * x))^2;

% Setup for the optimization.
p = p0;
ob1_plot(x0, p, dims, 'quick', 0);

% Optimize!
for k = 1 : num_iters

        %
        % Solve sub-problem for field.
        %

    % A_x = A{1} * D_(A{2} * p) * A{3} - omega^2 * speye(N) + eta * D_(env(:)); 
    A_x = S.r' * (A{1} * D_(A{2} * p) * A{3} - omega^2 * speye(N));

    % Try (slightly) different phases.
    x = (A_x * S.x) \ -(A_x * x0);
    x = S.x * x + x0;
    fprintf('%d: %e, ', k, phys_res(x, p));


        %
        % Solve sub-problem for structure.
        %

    A_p = A{1} * D_(A{3} * x) * A{2};
    b_p = omega^2 * x;

    A_hat = A_p * S.p;
    b_hat = b_p - A_p * p0;
    
    cvx_quiet(true);
    cvx_begin
        variable p(size(A_hat,2))
        minimize norm(A_hat * p - b_hat)
        subject to
            p >= 1/12.25 - S.p' * p0
            p <= 1 - S.p' * p0
    cvx_end
    p = S.p * p + p0;

    fprintf('%e\n', phys_res(x, p));

        %
        % Plot results and update.
        %

    ob1_plot(x, p, dims, 'quick', k);
end

% Export epsilon.
epsilon = reshape(1./p, dims);
