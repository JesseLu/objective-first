function ob1_demo1(dims)

N = prod(dims);
omega = 0.15;

    %
    % Form matrices.
    %

[A, S] = ob1_matrices(dims, dims/2);


    % 
    % Get initial value of p.
    %

% Form a simple waveguide.
epsilon = ones(dims);
epsilon(:,(dims-10)/2:(dims+10)/2) = 12.25;

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

D_ = @(z) sparse(1:length(z), 1:length(z), z, length(z), length(z));

A_x = A{1} * D_(A{2} * p0) * A{3} - omega^2 * speye(N); 
x = (A_x * S.x) \ -(A_x * x0);
x = S.x * x + x0;

A_p = A{1} * D_(A{3} * x) * A{2};
b_p = omega^2 * x;
p = (A_p * S.p) \ (b_p - A_p * p0);
p = S.p * p + p0;
ob1_plot(x, p, dims, 'quick');
