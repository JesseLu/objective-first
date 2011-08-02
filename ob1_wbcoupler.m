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
%         The initial structure.
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



