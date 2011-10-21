function [Ex, Ey, Hz] = ob1_fdfd(spec, eps, t_pml, pad)
% 
% Description
%     Solve an FDFD (finite-difference, frequency-domain) problem.

sigma_pml = 1 / spec.omega; % Strength of PML.
exp_pml = 3.5; % Exponential spatial increase in pml strength.
dims = size(eps);
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

b = zeros(dims); % Input excitation, equivalent to magnetic current source.
% in_pos = max([t_pml+1, round(pad(1)/2)]); % Location of input excitation.
in_pos = 5;

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
b(in_pos+1, pad(3)+1:end-pad(4)) = spec.in.Hz;
b(in_pos, pad(3)+1:end-pad(4)) = -spec.in.Hz * exp(i * spec.in.beta);

b = b ./ eps_y; % Convert from field to current source.

% Normalization factor so that the input power is unity.
b = -i * 2 * spec.in.beta / (1 - exp(i * 2 * spec.in.beta)) *  b;

b = b(:); % Vectorize.


    %
    % Solve.
    %

Hz = A \ b; % This should be using sparse matrix factorization. 

E = 1/spec.omega * inv_eps * Hcurl * Hz; % Obtain E-fields.

% Reshape and extract all three fields.
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);

