function [x, x_update] = ob1_field_setup(omega, phi, phi2eps, in, out, ...
    initial_option, update_option)
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
eps = 1./(A_spread * phi2eps(phi));
eps = struct('x', reshape(eps(1:N), dims), 'y', reshape(eps(N+1:2*N), dims));


% Solve for the input and output waveguide modes.
x = ob1_priv_wgmode(omega, eps, in{1}, in{2}) + ...
    ob1_priv_wgmode(omega, eps, out{1}, out{2});


    %
    % Form the physics residual function.
    %

phys_res = @(x, phi) norm(A(phi2eps(phi)) * x)^2;


    %
    % Form the field update function based on gradient descent.
    %

% Create template to hold border values constant.
tp = zeros(dims);
tp(3:end-2,3:end-2) = 1;
template = repmat(tp(:), 1, 1);

% Form the update x function.
x_update = @(x, phi) my_update_x(A(phi2eps(phi)), template, x, update_option);

% Determine the initial value for x.
switch initial_option
    case 'border-vals' % Start with only the boundary values.
    case 'soft-solve' % Start with a full soft-solve of the field.
        x = my_update_x(A(phi2eps(phi)), template, x, 'optimal');
    otherwise
        error('Invalid option for initial field.');
end


function [x, res] = my_update_x(A, P, x, update_option)

switch update_option
    case 'gradient'
        r = A * x; % Residual.
        g = P .* (A' * r); % Gradient.
        h = A * g; % Used to calculate step.
        s = (g'*g) / (h'*h); % Optimal step size (may have numerical error).

        x = x - s * g; % Update x.
        res = norm(A*x)^2; % New residual.
    case 'optimal'
        % Create selection matrix for active elements of x.
        ind = find(P);
        m = length(ind);
        S = sparse(1:m, ind, ones(m, 1), m, length(x));

        x0 = ~P .* x; % Constant border of x.

        x = (A * S') \ (-A * x0);
        x = S' * x + x0;

        res = norm(A*x)^2; % New residual.
end
