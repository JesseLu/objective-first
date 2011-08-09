function [phi, residuals] = ob1_wgcoupler(omega, epsilon, active_box)
% [PHI, RESIDUALS] = OB1_WGCOUPLER(OMEGA, EPSILON, ACTIVE_BOX)
% 
% Description
%     Perform an objective-first design of a waveguide-to-waveguide coupler.
% 
%     OB1_WGCOUPLER attempts to design a waveguide-to-waveguide coupler by
%     modifying the structure within a centered box within the simulation grid.
% 
%     The input and output waveguide modes are determined as the fundamental 
%     modes on the left and the right of the grid respectively. To calculate
%     these modes only the structure parameters at the very left and right 
%     edges of the grid are used. 
% 
%     Once the input and output modes are determined, the electromagnetic 
%     fields and dielectric structure, are alternately updated to decrease the
%     physics residual.
% 
% Inputs
%     OMEGA: Positive scalar.
%         Unitless frequency.
% 
%     EPSILON: 2D array.
%         The initial structure, which must be composed of only two materials.
%        
%     ACTIVE_BOX: 2-element array.
%         ACTIVE_BOX(1) and ACTIVE_BOX(2) describe the width and height of the
%         box in which the structure is allowed to vary. The box is centered
%         within the grid.
% 
% Outputs
%     PHI: 2D array (level-set function).
%         Level-set function describing the optimized structure.
% 
%     RESIDUALS: Vector.
%         The value of the physics residual at every sub-iteration.
% 
% Notes
%     Requires both lset-opt and wave-tools packages, both of which can be
%     found at github.com/JesseLu.


% Add access to the lset-opt and the wave-tools packages.
path(path, '~/lset-opt/matlab');
path(path, '~/wave-tools/em_bval_2dte');
path(path, '~/wave-tools/helper');

dims = size(epsilon);


    % 
    % Structure setup including:
    %   * converting EPSILON to a valid phi (level-set function),
    %   * determining the limiting values of EPSILON,
    %   * forming the phi2eps conversion function, and
    %   * forming the update structure function.
    %

[phi, phi2eps, phi_update] = ob1_structure_setup(omega, epsilon, active_box);


    %
    % Field setup including:
    %   * determining the input and output waveguide modes, and
    %   * forming the update field function.
    %

% The input is the fundamental incoming waveguide mode on the left, and the
% outgoing mode is the fundamental outgoing mode on the right.
[x, x_update] = ob1_field_setup(omega, phi, phi2eps, ...
    {'x-', 'in', 1}, {'x+', 'out', 1});


    %
    % Perform the field and structure optimization.
    %

ob1_plot(x, phi2eps(phi), dims, 'quick');

for k = 1 : 1e3
    [x, x_res(k)] = x_update(x, phi);
    [phi, phi_res(k)] = phi_update(x, phi);
    % lso_plot(phi); pause(0.1);
    if (mod(k, 100) == 1)
        ob1_plot(x, phi2eps(phi), dims, 'quick');
    end
end
ob1_plot(x, phi2eps(phi), dims, 'full');
figure(3); semilogy([x_res; phi_res]');

% phi_update(x, phi)
