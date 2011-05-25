function [A, b, E2H] = ...
    setup_physics(dims, omega, t_pml, p2e, e2p)

N = prod(dims);

sigma = 1 / omega; % Strength of PML.


    %
    % Helper functions for building matrices.
    %

% Allows for a clean definition of curl.
global S D

% Stretched-coordinate PML absorbing layers.
scx = @(sx, sy) setup_stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma, 2.5);
scy = @(sx, sy) setup_stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma, 2.5);

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   scy(.5,.5)*-(S(0,1)-S(0,0)),   scx(.5,.5)*(S(1,0)-S(0,0))];  
Hcurl = [   scy(.5,0)*(S(0,0)-S(0,-1));      scx(0,.5)*-(S(0,0)-S(-1,0))]; 

    % 
    % Build the matrices that we will be using.
    %

% Primary physics matrix, electromagnetic wave equation.
A = @(eps) Hcurl * Ecurl - omega^2 * D([eps.x(:); eps.y(:)]);

% Source terms.
b = @(J, M) (i * omega * J(:) - Hcurl * M(:));

E2H = @(E) 1 / (i * omega) * Ecurl * E;
