function [epsilon] =  ob1_demo1(dims, num_iters)

N = prod(dims);
omega = 0.15;
path(path, genpath('~/lumos/cvx'));

    %
    % Form matrices.
    %

[A, S] = ob1_matrices(dims, dims/2);


    % 
    % Get initial value of p.
    %

% Form a simple waveguide.
epsilon = ones(dims);
epsilon(:,(dims(2)-10)/2:(dims(2)+10)/2) = 12.25;
epsilon((dims(1)-10)/2:(dims(1)+10)/2,:) = 12.25;
epsilon(dims(1)/3:2*dims(1)/3, :) = 12.25;

% p0 = 1 ./ (epsilon(:) - (S.p * S.p') * (epsilon(:) - 1));
p0 = 1 ./ epsilon(:);

    
    % 
    % Get initial value of x.
    %

eps = A{2} * epsilon(:);
eps = struct('x', reshape(eps(1:N), dims), 'y', reshape(eps(N+1:2*N), dims));

x0 = ob1_wgmode(omega, eps, 'x-', 'in') + ...
    ob1_wgmode(omega, eps, 'x+', 'out');


    %
    % Optimize.
    %

% Helper function for sparse diagonal matrix.
D_ = @(z) sparse(1:length(z), 1:length(z), z, length(z), length(z));

% Physics residual calculation.
phys_res = @(x, p) norm(S.r' * (A{1} * D_(A{2} * p) * A{3} * x - omega^2 * x))^2;

p = p0;
for k = 1 : num_iters

    % Solve sub-problem for field.
    A_x = A{1} * D_(A{2} * p) * A{3} - omega^2 * speye(N); 
    x = (A_x * S.x) \ -(A_x * x0);
    x = S.x * x + x0;
    fprintf('%d: %e, ', k, phys_res(x, p));

    % Solve sub-problem for structure.
    A_p = A{1} * D_(A{3} * x) * A{2};
    b_p = omega^2 * x;

    A_hat = A_p * S.p;
    b_hat = b_p - A_p * p0;
    % p = (A_p * S.p) \ (b_p - A_p * p0);
    
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
    ob1_plot(x, p, dims, 'quick');
end

epsilon = reshape(1./p, dims);
