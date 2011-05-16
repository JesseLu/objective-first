% RESONATOR_DEMO
% 
% Objective-first optimization of a 2D TE nanophotonic resonator.
help resonator_demo


    %
    % Some optimization parameters.
    %

dims = [60 40]; % Size of the grid.
N = prod(dims);

omega = 0.9; % Angular frequency of desired mode.


    % 
    % Build the matrices that we will be using.
    %

% Allows for a clean definition of curl.
s = @(sx, sy) mirror_shift(dims, -[sx sy]); % Mirror boundary conditions.

% Define the curl operators.
Ecurl = [-(s(0,0)-s(0,-1)), s(0,0)-s(-1,0)]; % Applied to E.
Hcurl = [s(0,1)-s(0,0); -(s(1,0)-s(0,0))]; % Applied to H. 

J = zeros(2*N, 1);
M = zeros(N, 1);
J(round(dims(1)/2 + dims(1)*dims(2)/2)) = i;
mu_inv = speye(N);
A = Hcurl * Ecurl - omega^2 * speye(2*N);
b = i * omega * J - Hcurl * mu_inv * M;
E = A \ b;
% [E, D] = eigs(A, 1, 'SM');
H = Ecurl * E;
plot_fields(dims, {'Ex', E(1:N)}, {'Ey', E(N+1:end)}, {'Hz', H});
