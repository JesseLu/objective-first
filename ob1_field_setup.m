function [x, x_update] = ob1_field_setup(omega, phi, phi2eps, in, out)
% [X, X_UPDATE] = OB1_FIELD_SETUP(OMEGA, PHI, PHI2EPS, IN, OUT)
% 
% Description
%     Determine the initial field, as well as the update field function. 
% 
% Inputs
%     OMEGA: Scalar (unitless frequency).
% 
%     PHI: 2D array (level-set function).
%         Describes the structure.
% 
%     PHI2EPS: Function handle.
%         Convert PHI to epsilon values.
% 
%     IN, OUT: 3-element cell.
%         Determine the edge ('x-', 'x+', 'y-', or 'y+'), directionality 
%         ('in' or 'out'), and mode order (1 = fundamental mode, 2 = second-
%         order mode, etc..) of the input and output waveguide modes.
% 
% Outputs
%     X: Vector (Initial value of field).
% 
%     X_UPDATE: Function handle.
%         X_UPDATE(X, PHI) returns X such that the physics residual is 
%         decreased.

dims = size(phi);
N = prod(dims);

    %
    % Determine input and output waveguide modes, based on structure.
    %

% Find epsilon.
[A, B, d, A_spread] = ob1_priv_physics(omega, dims);
eps = A_spread * phi2eps(phi);
eps = struct('x', reshape(eps(1:N), dims), 'y', reshape(eps(N+1:2*N), dims));


% Solve for the input and output waveguide modes.
x = ob1_priv_wgmode(omega, eps, in{1}, in{2}, 0) + ...
    ob1_priv_wgmode(omega, eps, out{1}, out{2}, 0);


    %
    % Form the physics residual function.
    %

phys_res = @(x, phi) norm(A(phi2eps(phi)) * x)^2;


    %
    % Form the field update function based on gradient descent.
    %

tp = ones(dims);
tp([1,dims(1)],:) = 0;
tp(:,[1,dims(2)]) = 0;
template = [tp(:); tp(:); tp(:)];
x_update = @(x, phi) my_update_x(A(phi2eps(phi)), template, x);
% fprintf('%e\n', phys_res(x, phi));


function [x, res] = my_update_x(A, P, x)

r = A * x; % Residual.
g = P .* (A' * r); % Gradient.
h = A * g; % Used to calculate step.
s = (g'*g) / (h'*h); % supposed optimal step size.

x = x - s * g;

res = norm(A*x)^2;
