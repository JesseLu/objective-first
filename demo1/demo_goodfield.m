% DEMO_GOODFIELD
% 
% The objective of this demo is to determine a good objective to optimize on.
help demo_goodfield


    %
    % Some optimization parameters.
    %

dims = [160 90]; % Size of the grid.
N = prod(dims);

eps_lo = 1.0; % Relative permittivity of air.
eps_hi = 12.25; % Relative permittivity of silicon.

omega = 0.2; % Angular frequency of desired mode.

pml_thick = 10; % Thickness of pml.

% Set the device boundary.
dev.dims = [100 60];
dev.offset = round((dims - dev.dims)/2);


    %
    % Helper function for determining derivative matrices.
    % Also, helper global variables for prettier argument passing.
    %

global S D MY_DIMS MY_DEVICE

% Shortcut to form a derivative matrix.
S = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D = @(x) spdiags(x(:), 0, numel(x), numel(x));

MY_DIMS = dims;
MY_DEVICE = dev;


    %
    % Form the initial structure.
    %

lset_grid(dims);
phi = lset_box([-80 -15], [1000 10]);
% phi = lset_union(phi, lset_box([80 15], [100 10]));
phi = lset_complement(phi);

lset_plot(phi); pause % Use to visualize the initial structure.
% Initialize phi, and create conversion functions.
[phi, phi2p, phi2e, p2e, e2p, phi_smooth] = ...
    setup_levelset(phi, eps_lo, eps_hi, 1e-3);


    %
    % Find the input and output modes.
    %

input = mode_solve(mode_cutout(phi2e(phi), '-x'), omega, '-x');
mode_insert(


    % 
    % Source terms.
    %

Jx = zeros(dims);
Jy = zeros(dims);
Jy(ceil(dims(1)/2), ceil(dims(2)/2)) = i;
J = [Jx, Jy];

M = zeros(dims);


    %
    % Get matrices.
    %

[A, b] = setup_physics(dims, omega, pml_thick, p2e, e2p);
b = b(J, M);
% This defines the design objective.
C = -i * J(:);
d = 1e3;


    % 
    % Stuff for optimizing the structure.
    %

% Physics residual.
phys_res = @(x, p) 0.5 * norm(A(p2e(p)) * x)^2;

% Gradient of the physics residual with respect to p.
grad_res = @(x, p) e2p((-omega^2 * D(x))' * (A(p2e(p)) * x)); 


    %
    % Plot field.
    %
x = field_update(A(p2e(ones(dims)))+0.1*speye(2*N), b, C, d);
plot_fields(dims, {'Ex', x(1:N)}, {'Ey', x(N+1:end)});
