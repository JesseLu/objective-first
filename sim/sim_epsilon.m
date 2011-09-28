function [Ex, Ey, Hz] = sim_epsilon(omega, epsilon, dir)

dims = size(epsilon.x);
N = prod(dims);

pml_thick = 10; % Thickness of pml.

% Set the device boundary.
dev.dims = round(dims/2);
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
    % Find the input and output modes.
    %

input = mode_solve(mode_cutout(epsilon, dir), omega, dir);
[Jx, Jy, M] = mode_insert(input, dir);
J = [Jx, Jy];


    %
    % Get the physics matrices and solve.
    %

% Obtain physics matrix.
[A, b, E2H] = setup_physics(dims, omega, pml_thick);

% Solve.
x = A(epsilon) \ b(J,M);

% Back-out field components.
Ex = reshape(x(1:N), dims);
Ey = reshape(x(N+1:end), dims);
Hz = reshape(E2H(x), dims);


    %
    % Plot results.
    %

% Plot the structure.
eps = epsilon;
figure(1); plot_fields(dims, {'\epsilon_x', eps.x}, {'\epsilon_y', eps.y});

% Plot the fields.
figure(2); plot_fields(dims, ...
    {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
    {'Im(Ex)', imag(Ex)}, {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
    {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});
