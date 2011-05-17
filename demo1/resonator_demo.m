% RESONATOR_DEMO
% 
% Objective-first optimization of a 2D TE nanophotonic resonator.
help resonator_demo


    %
    % Some optimization parameters.
    %

dims = [160 90]; % Size of the grid.
N = prod(dims);
omega = 0.2; % Angular frequency of desired mode.
eps_lo = 1.0; % Relative permittivity of air.
eps_hi = 12.25; % Relative permittivity of silicon.

    %
    % Construct the initial level set for the structure.
    %

lset_grid(dims);
phi = lset_complement(lset_circle([0 0], 10));
phi = signed_distance(phi, 1e-3); % Make phi roughly a signed distance function.
lset_plot(phi); % Use to visualize the initial structure.

% Helper function to convert from phi to p ("digitize phi").
p = @(phi) ((phi .* (phi > -1 & phi < 1) + (phi >= 1) - (phi <= -1)) * ... 
    (eps_hi-eps_lo)/2) + (eps_hi+eps_lo)/2;


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

[A, b, phys_res, grad_res] = setup_physics(dims, omega);
b = b(J, M);

% This defines the design objective.
C = -i * J(:);
d = 1e3;


    %
    % Iterate through the optimization process.
    %

x = field_update(A(p(phi)), b, C, d);
phys_res(x, p(phi))
V = derivatives(grad_res(x, p(phi)));
phi = update_interface(phi, V, 0); % Move the interface.
phys_res(x, p(phi))

plot_fields(dims, {'p', p(phi)}, {'Ex', x(1:N)}, {'Ey', x(N+1:end)});
