function demo2()
% DEMO2
%
% Update the structure using the c-go package.
help demo2

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
N = prod(dims);


    %
    % Make the initial structure.
    %

lset_grid(dims);
phi = lset_box([0 0], [1000 10]);
%  eta = lset_box([0 0], [40 40]);
%  phi = lset_intersect(phi, lset_complement(eta));
phi = lset_complement(phi);

% Initialize phi, and create conversion functions.
[phi2e, phi2eps, phi_smooth] = setup_levelset(phi, 1.0, 12.25, 1e-3);


    %
    % Obtain the field.
    %

% Get the matrix to solve for the field.
[A] = my_field_physics(omega, phi2eps(phi)); 

% Initial values.
[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));

% Obtain the matrices for the boundary-value problem.
[Ahat, bhat, add_border] = ...
    setup_border_insert(A, [Ex(:); Ey(:); Hz(:)]);

% Solve the boundary-value problem.
xhat = Ahat \ -bhat;

% Obtain the full field (re-insert the field values at the boundary).
x = add_border(xhat);


    %
    % Setup the constrained gradient optimization.
    %

% Objective function and its gradient.
[f, g] = my_structure_physics(omega, x);

% Setup constraints
tp = ones(dims);
tp([1,dims(1)],:) = 0;
tp(:,[1,dims(2)]) = 0;
tp = [tp(:); tp(:); tp(:)];

c = @(v, dv, s) struct('x', v.x - s * (tp .* dv.x));



    %
    % Optimize using the c-go package.
    %

[v, fval, ss_hist] = opt(f, g, c, v, 1e3);


    %
    % Plot results.
    %

Ex = reshape(v.x(1:N), dims);
Ey = reshape(v.x(N+1:2*N), dims);
Hz = reshape(v.x(2*N+1:3*N), dims);

figure(1); plot_fields(dims, ...
    {'Re(Ex)', real(Ex)}, {'Re(Ey)', real(Ey)}, {'Re(Hz)', real(Hz)}, ...
    {'Im(Ex)', imag(Ex)}, {'Im(Ey)', imag(Ey)}, {'Im(Hz)', imag(Hz)}, ...
    {'|Ex|', abs(Ex)}, {'|Ey|', abs(Ey)}, {'|Hz|', abs(Hz)});

figure(2); cgo_visualize(fval, ss_hist);



function [A] = my_field_physics(omega, eps)


    %
    % Helper functions for building matrices.
    %

global S_ DIMS_ D_
N = prod(DIMS_); 

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   -(S_(0,1)-S_(0,0)),  (S_(1,0)-S_(0,0))];  
Hcurl = [   (S_(0,0)-S_(0,-1)); -(S_(0,0)-S_(-1,0))]; 

e = [eps.x(:); eps.y(:)];

A = [Ecurl, -i*omega*speye(N); i*omega*D_(e), Hcurl];

function [f, g] = my_structure_physics(omega, x)
