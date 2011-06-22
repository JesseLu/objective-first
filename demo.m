function demo(dims, field_or_struct, cgo_iters)
% DEMO(DIMS, FIELD_OR_STRUCT, CGO_ITERS)
%
% Description
%     Independently test the optimization routines for the field and 
%     structure variables. Results are verified against "known" solutions.
% 
%     Basically, eye-ball it to see if the result from matrix-inversion makes
%     sense, then compare that (especially the value of the physics residual)
%     to the result given by gradient optimization.

help demo

path(path, '~/c-go'); % Make sure we have access to c-go.
path(path, '~/level-set'); % Make sure we have access to level-set.

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
eta = lset_box([0 0], dims/2);
% phi = lset_intersect(phi, lset_complement(eta));
phi = lset_complement(phi);

% Initialize phi, and create conversion functions.
[phi2p, phi2eps, phi_smooth] = setup_levelset(phi, 1.0, 12.25, 1e-3);


    %
    % Setup the constrained gradient optimization.
    %

% Objective function and its gradient.
[f, g] = em_physics2(omega); 

% Setup constraints
tp = ones(dims);
tp([1,dims(1)],:) = 0;
tp(:,[1,dims(2)]) = 0;
tp = [tp(:); tp(:); tp(:)];

% % This constraint function allows both variables to change.
% c = @(v, dv, s) struct( 'x', v.x - s * (tp .* dv.x), ...
%                         'p', v.p - s * ((eta(:) < 0) .* dv.p));

switch field_or_struct
    case 'field'
        c = @(v, dv, s) struct( 'x', v.x - s * (tp .* dv.x), ...
                                'p', v.p);
    case 'struct'
        c = @(v, dv, s) struct( 'x', v.x, ...
                                'p', v.p - s * ((eta(:) < 0) .* dv.p)); 
    otherwise
        error('Invalid parameter for FIELD_OR_STRUCT.');
end


% Initial values.
[Ex, Ey, Hz] = setup_border_vals({'x-', 'x+'}, omega, phi2eps(phi));
v.x = [Ex(:); Ey(:); Hz(:)];
% randn('state', 1);
% v.x = randn(size(v.x));
v.p = phi2p(phi);
v.p = v.p(:);
% v.p = randn(N, 1);


    %
    % Optimize by directly solving the matrix equation.
    %

        [A, b, reinsert] = em_physics1('field', omega, tp, v.x);
        v.x = reinsert(A(v.p) \ b(v.p));
switch(field_or_struct)
    case 'struct'
        [A, b, reinsert] = em_physics1('struct', omega, eta<0, v.p);
        v0.p = reinsert(A(v.x) \ b(v.x));
end

    %
    % Optimize using the c-go package.
    %

[v, fval, ss_hist] = opt(f, g, c, v, 1e2);


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

figure(2); plot_fields(dims, {'p', v0.p});

% figure(3); cgo_visualize(fval, ss_hist);


