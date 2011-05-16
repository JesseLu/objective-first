% RESONATOR_DEMO
% 
% Objective-first optimization of a 2D TE nanophotonic resonator.
help resonator_demo


    %
    % Some optimization parameters.
    %

dims = [60 40]; % Size of the grid.
N = prod(dims);

omega = 0.3; % Angular frequency of desired mode.

    % 
    % Build the matrices that we will be using.
    %

s = @(sx, sy) shift_grid(dims, -[sx sy]); % Allows for clean definition of curl.
curl = [s(0,0)-s(0,-1), s(0,0)-s(-1,0)]; % Here, dH/dt = curl E.

A = curl * curl' - omega^2 * speye(N);

% Test out by finding eigenmodes.
[H, D] = eigs(A, 1, 'SM');
E = curl' * H;
plot_fields(dims, {'Ex', E(1:N)}, {'Ey', E(N+1:end)}, {'Hz', H});
