% This is an example of an objective-first design of a nanophotonic waveguide 
% coupler.
help example

% Create the initial structure, a silicon waveguide in air.
eps0 = ones([40 40]);
eps0(:,10:30) = 12.25;

% Fill everything, except for a 2-cell boundary layer, with epsilon = 9.
% The 2-cell boundary layer is left unchanged because the boundary conditions
% use these cells.
% Feel free to comment out this line, you will still get a working design!
% eps0(3:end-2,3:end-2) = 9.0;

% Setup the design. Specifically, we want
% *   a design frequency of 0.15,
% *   to limit epsilon from 1 to 12.25, and
% *   to use the fundamental mode as an input on the left and
%     the second-order mode as the output mode on the right.
spec = setup(0.15, eps0, [1 12.25], [1 2]);
return
% Optimize for 40 iterations, or until the gradient norm drops below 1e-4.
eps = solve(spec, 40, 1e-6, 'alt_dir_mod');

% Simulate the result. Make the simulation space 100x100.
simulate(spec, eps, [160 120]);
