% RESONATOR_DEMO
% 
% Objective-first optimization of a 2D TE nanophotonic resonator.
help resonator_demo


    %
    % Some optimization parameters.
    %

dims = 2*[160 90]; % Size of the grid.
N = prod(dims);

omega = 0.5; % Angular frequency of desired mode.
sigma = 1e0/omega; % Strength of PML.
t_pml = 10; % Thickness of PML layer.

    % 
    % Build the matrices that we will be using.
    %

% Allows for a clean definition of curl.
s = @(sx, sy) mirror_shift(dims, -[sx sy]); % Mirror boundary conditions.
scx = @(sx, sy) stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma, 2.5);
scy = @(sx, sy) stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma, 2.5);

% Define the curl operators as applied to E and H, respectively.
Ecurl = [   scy(.5,.5)*-(s(0,0)-s(0,-1)),   scx(.5,.5)*(s(0,0)-s(-1,0))];  
Hcurl = [   scy(.5,0)*(s(0,1)-s(0,0));      scx(0,.5)*-(s(1,0)-s(0,0))]; 

J = zeros(2*N, 1);
M = zeros(N, 1);
J(prod(dims) + round(dims(1)/2 + dims(1)*dims(2)/2)) = i;
mu_inv = speye(N);
A = Hcurl * Ecurl - omega^2 * speye(2*N);
b = i * omega * J - Hcurl * mu_inv * M;
E = A \ b;
% [E, D] = eigs(A, 1, 'SM');
H = Ecurl * E;
plot_fields(dims, {'Ex', E(1:N)}, {'Ey', E(N+1:end)}, {'Hz', H});
