function [eff] = simulate(spec, eps, dims)
% EFF = SIMULATE(SPEC, EPS, DIMS)
% 
% Description
%     Simulates the design and determines the performance.
% 
% Inputs
%     SPEC: Structure.
%         This is the specification of the problem, as obtained from SETUP().
% 
%     EPS: 2-d array.
%         The permittivity values of the design to be simulated.
% 
%     DIMS: 2-element vector.
%         The size of the simulation in the x- and y- directions respectively.
%         These values should be considerable larger than size(EPS).
% 
% Outputs
%     EFF: Non-negative scalar.
%         The coupling efficiency where 1.0 is 100% conversion efficiency.
% 
% Examples

% Hard-coded parameters.
t_pml = 10; % Thickness of PML.
sigma_pml = 1 / spec.omega; % Strength of PML.
exp_pml = 3.5; % Exponential spatial increase in pml strength.

    
    % 
    % Determine epsilon for the simulation.
    %

% Number of extra cells to add to eps.
pad = [ floor((dims(1) - size(eps, 1))/2), ...
        ceil((dims(1) - size(eps, 1))/2), ...
        floor((dims(2) - size(eps, 2))/2), ...
        ceil((dims(2) - size(eps, 2))/2)];

% Expand eps to the full simulation size.
eps = cat(1, repmat(eps(1,:), pad(1), 1), eps, repmat(eps(end,:), pad(2), 1));
eps = cat(2, repmat(eps(:,1), 1, pad(3)), eps, repmat(eps(:,end), 1, pad(4)));

[eps_x, eps_y] = ob1_interp_eps(eps); % Get x and y components of eps.


    %
    % Build the simulation matrix.
    %

% Shortcut to form a derivative matrix.
S = @(sx, sy) ob1_shift_matrix(dims, -[sx sy]);

% Helper function to create stretched-coordinate PML absorbing layers.
scx = @(sx, sy) ob1_stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma_pml, exp_pml);
scy = @(sx, sy) ob1_stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma_pml, exp_pml);

% Define the curl operators as applied to E and H, respectively.
Ecurl = [scy(.5,.5)*-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
Hcurl = [scy(.5,0)*(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 

% Diagonal matrix for 1/epsilon.
inv_eps = spdiags([eps_x(:).^-1; eps_y(:).^-1], 0, 2*prod(dims), 2*prod(dims));

% This is the matrix that we will solve.
A = Ecurl * inv_eps * Hcurl - spec.omega^2 * speye(prod(dims));

    
    %
    % Determine the input excitation.
    %

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
b = zeros(dims);
b(round(dims(1)/2), pad(3)+1:end-pad(4)) = spec.in.Hz;
% b(round(pad(1)/2), pad(3)+1:end-pad(4)) = spec.in.Hz;
% b(round(pad(1)/2)-1, pad(3)+1:end-pad(4)) = -spec.in.Hz * ...
%                                             exp(i * spec.in.beta);

% Normalization factor so that the input power is unity.
b = spec.omega / (1 - exp(-i * 2 * spec.in.beta)) / 2 * b;


    %
    % Solve.
    %


Hz = A \ b(:);
Hz = reshape(Hz, dims);

    %
    % Calculate efficiency.
    %
        

max(abs(Hz(:)))

subplot 121; imagesc(abs(Hz)')
subplot 122; plot(abs(Hz(:,48)), '.-')
% subplot 121; plot([eps_x(dims(2)/2,:);eps_y(dims(2)/2,:)]', '.-');
