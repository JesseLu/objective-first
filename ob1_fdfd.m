function [Ex, Ey, Hz] = ob1_fdfd(omega, eps, in)
% 
% Description
%     Solve a FDFD (finite-difference, frequency-domain) problem, using the
%     input mode as the source term.
% 
%     OB1_FDFD actually expands the structure in order to add absorbing
%     boundary layers (pml), simulates the structure by sourcing it from
%     within the pml, and then truncates the solution so the user does not
%     see the absorbing pml region.
% 
%     Note that a time dependence of exp(-i * omega * t) is assumed for all
%     fields.


% Hard-coded parameters for the pml.
t_pml = 10; % Thickness of pml.
sigma_pml = 1 / omega; % Strength of pml.
exp_pml = 3.5; % Exponential spatial increase in pml strength.

% Expand eps to include room for the pml padding
eps = ob1_pad_eps(eps, t_pml * [1 1 1 1]);
dims = size(eps); % New size of the structure.
[eps_x, eps_y] = ob1_interp_eps(eps); % Obtain x- and y- components of eps.

    
    %
    % Form the system matrix which will be solved.
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
A = Ecurl * inv_eps * Hcurl - omega^2 * speye(prod(dims));

    
    %
    % Determine the input excitation.
    %

b = zeros(dims); % Input excitation, equivalent to magnetic current source.
in_pos = 2; % Cannot be 1, because eps interpolation wreaks havoc at border.

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
% b(in_pos+1, pad(3)+1:end-pad(4)) = in.Hz;
b_pad = [floor((dims(2) - length(in.Hz))/2), ...
        ceil((dims(2) - length(in.Hz))/2)];
b(in_pos, b_pad(1)+1:end-b_pad(2)) = -in.Hz * exp(i * in.beta);

b = b ./ eps_y; % Convert from field to current source.

b = b(:); % Vectorize.


    %
    % Solve.
    %

Hz = A \ b; % This should be using sparse matrix factorization. 

E = (i/omega) * inv_eps * Hcurl * Hz; % Obtain E-fields.

% Reshape and extract all three fields.
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);

% Cut off the pml parts.
Ex = Ex(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
Ey = Ey(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
Hz = Hz(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
