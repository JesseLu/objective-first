% DEMO1
%
% Update the field using the c-go package.
help demo1

path(path, '~/c-go'); % Make sure we have access to c-go.
path(path, '~/level-set'); % Make sure we have access to level-set.



dims = [80 80]; % Size of the grid.
omega = 0.15; % Angular frequency of desired mode.

    %
    % Helper function for determining derivative matrices.
    % Also, helper global variables for prettier argument passing.
    %

global S_ D_ DIMS_ 

% Shortcut to form a derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D_ = @(x) spdiags(x(:), 0, numel(x), numel(x));

DIMS_ = dims;


lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
phi = lset_complement(phi);

% Initialize phi, and create conversion functions.
[phi, phi2p, phi2e, phi2eps, p2e, e2p, p2eps, eps2p, ...
    A_spread, A_gather, phi_smooth] = ...
    setup_levelset(phi, 1.0, 12.25, 1e-3);

eps = phi2eps(phi);

[f, g] = em_physics(omega, phi2eps(phi));
[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));
v.E = [Ex(:), Ey(:)];
v.H = Hz(:);

% Setup constraints

figure(1); plot_fields(dims, ...
    {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
    {'Im(Ex)', imag(Ex)}, {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
    {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});

return
N = prod(dims);

eps_lo = 1.0; % Relative permittivity of air.
eps_hi = 12.25; % Relative permittivity of silicon.



[A0, B, d] = setup_physics(omega, A_spread);
return
    %
    % Helper function for determining derivative matrices.
    % Also, helper global variables for prettier argument passing.
    %

global S_ D_ DIMS_ 

% Shortcut to form a derivative matrix.
S_ = @(sx, sy) shift_mirror(dims, -[sx sy]); % Mirror boundary conditions.

% Shortcut to make a sparse diagonal matrix.
D_ = @(x) spdiags(x(:), 0, numel(x), numel(x));

DIMS_ = dims;


    %
    % Form the initial structure.
    %

lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
% phi = lset_intersect(phi, lset_complement(lset_box([0 0], [3 1000]))); % Gap.
phi = lset_union(phi, lset_box([0 0], [10 1000])); % Cross-beam.
phi = lset_complement(phi);

% Initialize phi, and create conversion functions.
[phi, phi2p, phi2e, phi2eps, p2e, e2p, p2eps, eps2p, ...
    A_spread, A_gather, phi_smooth] = ...
    setup_levelset(phi, eps_lo, eps_hi, 1e-3);

% lset_plot(phi); pause % Use to visualize the initial structure.


    %
    % Find the input and output modes.
    %

[Ex, Ey, Hz] = setup_border_vals({'x-', 'y-'}, omega, phi2eps(phi));

    %
    % Mark the region where we will allow the structure to change.
    %

eta = lset_box([0 0], [40 40]);

p = phi2p(phi);
for k = 1 : 100

        %
        % Get the physics matrices and solve.
        %

    % Obtain physics matrix.
    [A0, B, d] = setup_physics(omega, A_spread);
    A = @(phi) A0(phi2p(phi));

    % Obtain the matrices for the boundary-value problem.
    [Ahat, bhat, add_border, rm_border] = ...
        setup_border_insert(A0(p), [Ex(:); Ey(:); Hz(:)]);

    % Solve the boundary-value problem.
    xhat = Ahat \ -bhat;

    % Obtain the full field (re-insert the field values at the boundary).
    x = add_border(xhat);

    % Calculate the values for the physics residual.
    phys_res = @(p, x) norm(rm_border(A0(p)*x));
    % norm(Ahat*xhat + bhat) % Alternate definition.

    fprintf('%d: %e', k, phys_res(p, x));

    p = solve_p(B(x), d(x), phi2p(phi), eta);
    fprintf(' -> %e\n', phys_res(p, x));
    Hz = reshape(x(2*N+1:end), dims);
    plot_fields(dims, {'p', p}, {'|Hz|', abs(Hz)});
    drawnow
end

    %
    % Plot results.
    %

% Back-out field components.
Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:2*N), dims);
Hz = reshape(x(2*N+1:end), dims);

% Plot the structure.
eps = p2eps(p);
figure(1); plot_fields(dims, {'\epsilon_x', eps.x}, {'\epsilon_y', eps.y});

% Plot the fields.
figure(2); plot_fields(dims, ...
    {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
    {'Im(Ex)', imag(Ex)}, {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
    {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});

% % Plot cross-section, and check power flow.
% ey = Ey(:, dims(2)/2);
% hz = Hz(:, dims(2)/2);
% figure(3); plot([real(ey), real(hz), real(conj(hz).*ey)]);
% figure(3); plot_fields(dims, {'power', real(Ex.*conj(Hz))});
figure(3); plot_fields(dims, {'p', p});
drawnow

% % Plot the residual.
% y = real(A(phi)*x);
% y = rm_border(y);
% M = prod(dims-2);
% figure(4); plot_fields(dims-2, {'residual_x', y(1:M)}, ...
%     {'res_y', y(M+1:2*M)}, {'res_z', y(2*M+1:3*M)});
