function [A, b, phys_res, grad_res] = setup_physics(dims, omega)

N = prod(dims);

sigma = 1e0/omega; % Strength of PML.
t_pml = 10; % Thickness of PML layer.


    %
    % Helper functions for building matrices.
    %

% Allows for a clean definition of curl.
s = @(sx, sy) mirror_shift(dims, -[sx sy]); % Mirror boundary conditions.

% Stretched-coordinate PML absorbing layers.
scx = @(sx, sy) stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma, 2.5);
scy = @(sx, sy) stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma, 2.5);

% Shortcut to make a sparse diagonal matrix.
my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x));

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   scy(.5,.5)*-(s(0,0)-s(0,-1)),   scx(.5,.5)*(s(0,0)-s(-1,0))];  
Hcurl = [   scy(.5,0)*(s(0,1)-s(0,0));      scx(0,.5)*-(s(1,0)-s(0,0))]; 

% Epsilon is defined on Ez, interpolate to get epsilon at Ex and Ey.
A_spread = 0.5 * [s(0,0)+s(1,0); s(0,0)+s(0,1)];
spread_eps = @(p) my_diag(A_spread * p(:));

% Average Ex and Ey, to get back p (Ez).
A_gather = 1/4 * [s(0,0)+s(-1,0), s(0,0)+s(0,-1)];
gather_eps = @(eps) reshape(A_gather * eps(:), dims);


    % 
    % Build the matrices that we will be using.
    %

% Primary physics matrix, electromagnetic wave equation.
A = @(p) Hcurl * Ecurl - omega^2 * spread_eps(p);

% Source terms.
b = @(J, M) i * omega * J(:) - Hcurl * M(:);

% Physics residual.
phys_res = @(x, p) 0.5 * norm(A(p) * x)^2;

% Gradient of the physics residual with respect to p.
grad_res = @(x, p) gather_eps((-omega^2 * my_diag(x))' * (A(p) * x)); 


