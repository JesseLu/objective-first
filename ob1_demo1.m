function [epsilon] =  ob1_demo1(omega, epsilon, eps_lims, out_dir, num_iters, options)

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
p0 = (epsilon(:) - (S.p * S.p') * (epsilon(:) - eps_lims(2)));

    
    % 
    % Get initial value of x.
    %

eps = A.spread * epsilon(:);
eps = struct('x', reshape(eps(1:N), dims), 'y', reshape(eps(N+1:2*N), dims));

mode = {ob1_wgmode(omega, eps, 'x-', 'in'), ...
    ob1_wgmode(omega, eps, out_dir, 'out')};
x0 = mode{1} + mode{2};


    %
    % Optimize.
    %

% Helper function for sparse diagonal matrix.
D_ = @(z) sparse(1:length(z), 1:length(z), z, length(z), length(z));

% Physics residual calculation.
phys_res = @(x, p) ...
    norm(S.r' * (A.hcurl * A.ecurl * x - omega^2 * D_(A.spread * p) * x))^2;

% Setup for the optimization.
x = x0;
p = p0;
eta = options(1);
deta = options(2);

ob1_plot(x, p, dims, 'quick');

% Optimize!
for k = 1 : num_iters
        %
        % Solve sub-problem for field.
        %

    % A_x = A{1} * D_(A{2} * p) * A{3} - omega^2 * speye(N) + eta * D_(env(:)); 
    % A_x = S.r' * (A{1} * D_(A{2} * p) * A{3} - omega^2 * speye(N));
    A_x = S.r' * (A.hcurl * A.ecurl - omega^2 * D_(A.spread * p));
    A_div = S.d' * A.div * D_(A.spread * p);
    % A_x = [A_x, A.div'; A.div, sparse(N, N)];

    % x = (A_x * S.x) \ -(A_x * x0);
    x = ob1_priv_lseq(A_x * S.x, -(A_x * x0), ...
        A_div * S.x, -A_div * x0, eta);
    x = S.x * x + x0;
    fprintf('%d: (%e) %e, ', k, norm(A_div * x), phys_res(x, p));


        %
        % Solve sub-problem for structure.
        %

    A_p = omega^2 * D_(x) * A.spread;
    b_p = A.hcurl * A.ecurl * x;

    A_hat = A_p * S.p;
    b_hat = b_p - A_p * p0;
    
    cvx_quiet(true);
    cvx_begin
        variable p(size(A_hat,2))
        minimize norm(A_hat * p - b_hat)
        subject to
            p <= 12.25 - S.p' * p0
            p >= 1 - S.p' * p0
    cvx_end
    p = S.p * p + p0;

    fprintf('%e\n', phys_res(x, p));

        %
        % Plot results.
        %

    ob1_plot(x, p, dims, 'quick');
    eta = eta * deta;

end

% Export epsilon.
epsilon = reshape(1./p, dims);
