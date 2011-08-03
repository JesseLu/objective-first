function [phi, phi2eps, phi_update] = ...
    ob1_structure_setup(omega, epsilon, active_box)
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

       
    %
    % Convert epsilon to a valid phi (level-set function).
    %

eps_lims = [min(epsilon(:)), max(epsilon(:))]; % Find min and max epsilon.
phi = lso_regularize(epsilon - mean(eps_lims));


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

phys_res = @(x, phi) norm(B(x) * phi2eps(phi) - d(x))^2;
phi_update = phys_res;
