% DEMO2
%
% Update both structure and field in an unbounded way to satisfy physics.
help demo2


    %
    % Some optimization parameters.
    %

dims = [80 80]; % Size of the grid.
N = prod(dims);

eps_lo = 1.0; % Relative permittivity of air.
eps_hi = 12.25; % Relative permittivity of silicon.

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


    %
    % Form the initial structure.
    %

lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
phi = lset_intersect(phi, lset_complement(lset_box([0 0], [10 1000]))); % Form cross-beam.
phi = lset_complement(phi);

% Initialize phi, and create conversion functions.
[phi, phi2p, phi2e, phi2eps, p2e, e2p, p2eps, eps2p, ...
    A_spread, A_gather, phi_smooth] = ...
    setup_levelset(phi, eps_lo, eps_hi, 1e-3);

% lset_plot(phi); pause % Use to visualize the initial structure.


    %
    % Find the input and output modes.
    %

[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));


    %
    % Get the physics matrices and solve.
    %

% Obtain physics matrix.
[A0, B, d] = setup_physics(omega, p2e);
A = @(phi) A0(phi2p(phi));

% Obtain the matrices for the boundary-value problem.
[Ahat, bhat, add_border, rm_border] = ...
    setup_border_insert(A(phi), [Ex(:); Ey(:); Hz(:)]);

% Solve the boundary-value problem.
xhat = Ahat \ -bhat;

% Obtain the full field (re-insert the field values at the boundary).
x = add_border(xhat);

% Back-out field components.
Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:2*N), dims);
Hz = reshape(x(2*N+1:end), dims);

% Calculate the values for the physics residual.
phys_res = @(phi, x) norm(rm_border(A(phi)*x));
fprintf('Physics residual: %e\n', phys_res(phi, x));
% norm(Ahat*xhat + bhat) % Alternate definition.


    %
    % Mark the region where we will allow the structure to change.
    %

eta = lset_box([0 0], [40 40]);


    % 
    % Compute the partial derivative of the physics residual relative to p.
    %

p = [real(B(x)); imag(B(x))] * A_spread \ [real(d(x)); imag(d(x))];
% p0 = phi2p(phi);
% p = p0 .* (eta >= 0) + ...
%     e2p(real((conj(B(x)) .* d(x))) ./ abs(B(x)).^2) .* (eta < 0);
ind = find(isnan(p));
p(ind) = p0(ind);
norm(rm_border(A0(p) * x))
figure(3); plot_fields(dims, {'p', p});
% res = [];
% p = phi2p(phi);
% for k = 1 : 100
%     dp = (eta < 0) .* real(e2p(B(x)' * (B(x) * phi2e(phi) - d(x))));
%     p = p - 100 * dp;
%     % phi = update_interface(phi, [], real(dp), 0, 1e-2);
%     res(k) = norm(rm_border(A0(p) * x));
%     fprintf('.');
% end
% figure(3); subplot 111; plot(res, '.-')
% min(res)/max(res)

% dp = (eta < 0) .* e2p(B(x)' * (B(x) * phi2e(phi) - d(x)));
% c = -1 : 1e-2 : 1;
% for k = 1 : length(c) 
%     phi1 = update_interface(phi, [], real(dp), 0, c(k));
%     res(k) = phys_res(phi1, x);
% end
% figure(3); subplot 111; plot(c, res, '.-')


    %
    % Plot results.
    %

% Plot the structure.
eps = phi2eps(phi);
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

% % Plot the residual.
% y = real(A(phi)*x);
% y = rm_border(y);
% M = prod(dims-2);
% figure(4); plot_fields(dims-2, {'residual_x', y(1:M)}, ...
%     {'res_y', y(M+1:2*M)}, {'res_z', y(2*M+1:3*M)});
