function [x] = ob1_priv_wgmode(omega, eps, edge, in_out, phase)
% 
% Description
%     Solve for a particular waveguide mode at the edge of the grid.


global DIMS_
DIMS_ = size(eps.x);

% Find mode.
mode = mode_solve(mode_cutout(eps, edge), omega, edge);

% Insert single mode into empty 2D array.
[Ex, Ey, Hz] = mode_insert(mode, edge, in_out, phase);

x = [Ex(:); Ey(:); Hz(:)];

clear global DIMS_
