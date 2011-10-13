% This is an example of an objective-first design of a nanophotonic waveguide 
% coupler.
help example

% Create the initial structure, a silicon waveguide in air.
eps0 = ones([40 60]);
eps0(:,20:40) = 12.25;
eps0(3:end-2,3:end-2) = 9.0; 

% Setup the design. Specifically, we want
% *   a design frequency of 0.15,
% *   to limit epsilon from 1 to 12.25, and
% *   to use the fundamental mode as an input on the left and
%     the second-order mode as the output mode on the right.
spec = setup(0.15, eps0, [1 12.25], [1 2]);

% Optimize for 40 iterations, or until the gradient norm drops below 1e-4.
eps = solve(spec, 400, 1e-6);

% Simulate the result. Make the simulation space 100x100.
simulate(spec, eps, [200 100]);
