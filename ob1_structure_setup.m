function [phi, phi2eps, phi_update] = ob1_structure_setup(omega, epsilon, ...
    active_box, initial_option, update_option)
% [PHI, PHI2EPS, PHI_UPDATE] = OB1_STRUCTURE_SETUP(OMEGA, EPSILON, ACTIVE_BOX)
% 
% Description
%     Setup the structure for objective-first optimization. This involves:
%         * converting epsilon to a valid phi,
%         * forming the phi2eps conversion function, and
%         * forming the phi_update function.
% 
% Inputs
%     OMEGA: Scalar (unitless frequency).
% 
%     EPSILON: 2D array.
%         Initial values of permittivity.
% 
%     ACTIVE_BOX: 2-element array.
%         Determines the width and the height of the centered region where
%         the structure may be modified.
% 
% Outputs
%     PHI: 2D array (level-set function).
%         Initial level-set function describing the structure.
% 
%     PHI2EPS: Function handle.
%         PHI2EPS accepts a level-set function as input, and returns the values
%         of epsilon (permittivity) on the grid.
% 
%     PHI_UPDATE: Function handle.
%         PHI_UPDATE(X, PHI) returns an updated level-set function which 
%         decreases the physics residual of the system.

dims = size(epsilon);
N = numel(epsilon);

    % 
    % Active box isolation template.
    %

template = zeros(dims);
template(round((dims(1)-active_box(1))/2):round((dims(1)+active_box(1))/2), ...
    round((dims(2)-active_box(2))/2):round((dims(2)+active_box(2))/2)) = 1;

       
    %
    % Convert epsilon to a valid initial phi (level-set function).
    %

eps_lims = [min(epsilon(:)), max(epsilon(:))]; % Find min and max epsilon.

% Modify the initial epsilon according to various options.
switch initial_option
    case 'full-struct' % Keep the full structure (don't change it).
    case 'empty-box' % Empty active box.
        epsilon = template * eps_lims(1) + (~template) .* epsilon; 
    case 'cheese' % Put cheese in the active box.
        epsilon = template .* (mean(eps_lims) + lso_cheese(dims)) + (~template) .* epsilon; 
    otherwise
        error('Invalid option for initial structure.');
end

phi = lso_regularize(epsilon - mean(eps_lims)); % Convert to level-set function.


    %
    % Form phi to epsilon conversion function.
    %

[A, B, d, A_spread] = ob1_priv_physics(omega, size(phi));
phi2eps = @(phi) ...
    reshape(diff(eps_lims)/2 * lso_fracfill(phi) + mean(eps_lims), ...
    numel(phi), 1);

       
    %
    % Form update function for phi.
    %

phi_update = @(x, phi) my_phi_update(B(x), d(x), template, phi, phi2eps, ...
    update_option);



function [phi, res] = my_phi_update(B, d, P, phi, phi2eps, update_option)

phys_res = @(phi) norm(B * phi2eps(phi) - d)^2;
switch update_option
    case 'fixed' % Don't change phi.
        res = phys_res(phi);
    case 'unbounded' % Gradient-descent, without level-sets.
    case 'level-set' % Use level-sets.
        r = B * phi2eps(phi) - d; % Residual.
        g = real(P(:) .* (B' * r)); % Gradient.
        h = B * g; % Used to calculate step.
        s = (g'*g) / abs(h'*h); % Optimal step size (may have numerical error).

        dp = reshape(-s * g, size(phi));

        phi = lso_update(phi, 1e0 * dp, phys_res, 2.^[-40:10], P);
        phi = lso_quickreg(phi);
        res = phys_res(phi);
        fprintf('%e -> %e\n', norm(r)^2, res);
    otherwise
        error('Invalid option for structure update.');
end
