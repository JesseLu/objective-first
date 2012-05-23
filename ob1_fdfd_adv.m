function [Ex, Ey, Hz] = ob1_fdfd_adv(omega, eps, mu, input_exc, bc, t_pml)
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

dims = size(mu);
    
    %
    % Form the system matrix which will be solved.
    %

% Shortcut to form a derivative matrix.
S = @(sx, sy) ob1_shift_matrix(dims, -[sx sy]);

% Hard-coded parameters for the pml.
sigma_pml = 1 / omega; % Strength of pml.
exp_pml = 2.5; % Exponential spatial increase in pml strength.


% Helper function to create stretched-coordinate PML absorbing layers.
scx = @(sx, sy) ob1_stretched_coords(dims, [1 dims(1)+0.5], [sx, sy], ...
    'x', t_pml, sigma_pml, exp_pml);
scy = @(sx, sy) ob1_stretched_coords(dims, [1 dims(2)+0.5], [sx, sy], ...
    'y', t_pml, sigma_pml, exp_pml);


% Define the curl operators as applied to E and H, respectively.
if strcmp(bc, 'pml')
    Ecurl = [scy(.5,.5)*-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
    Hcurl = [scy(.5,0)*(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 
elseif strcmp(bc, 'per')
    Ecurl = [-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
    Hcurl = [(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 
end

% Diagonal matrix for 1/epsilon.
inv_eps = spdiags([eps{1}(:).^-1; eps{2}(:).^-1], 0, 2*prod(dims), 2*prod(dims));
mu = spdiags(mu(:), 0, prod(dims), prod(dims));

% This is the matrix that we will solve.
A = Ecurl * inv_eps * Hcurl - omega^2 * mu;

    
    %
    % Determine the input excitation.
    %

b = input_exc(:); % Vectorize.


    %
    % Solve.
    %

Hz = A \ b; % This should be using sparse matrix factorization. 

E = (i/omega) * inv_eps * Hcurl * Hz; % Obtain E-fields.

% Reshape and extract all three fields.
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);

