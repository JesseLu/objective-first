function [Ex, Ey, Hz, gradient] = ob1_fdfd_adj(omega, eps, in, g_x_partial, bc)
% 
% Description
%     Solve a transverse magnetic FDFD (finite-difference,
%     frequency-domain) problem, using the input mode as the source term.
%     Also, compute the gradient with respect to the permittivity (eps).
% 
%     OB1_FDFD actually expands the structure in order to add absorbing
%     boundary layers (pml), simulates the structure by sourcing it from
%     within the pml, and then truncates the solution so the user does not
%     see the absorbing pml region.
% 
%     Note that a time dependence of exp(-i * omega * t) is assumed for all
%     fields.
%
%     For objective function g computed on the field values, x, pass in the
%     function handle g_x_partial, which computes the partial deriviative
%     of the objective function with respect to the fields, x.  Note, given
%     that the objective function output is assumed to be real, but the
%     variable x is complex, just pass in the derivative with respect to x,
%     treating conj(x) as a separate variable.  Inside this function, we
%     take care of the resulting twice the real part of the result as is
%     necessary when adding the conj(x) derivative to the computation of
%     the derivative with respect to x.
%     

% Hard-coded parameters for the pml.
t_pml = 20; % Thickness of pml.
sigma_pml = 1 / omega; % Strength of pml.
exp_pml = 2.5; % Exponential spatial increase in pml strength.

% Get the unpadded size of epsilon
[height, width] = size(eps);

% Expand eps to include room for the pml padding
if strcmp(bc, 'pml') 
    eps = ob1_pad_eps(eps, t_pml * [1 1 1 1]);
elseif strcmp(bc, 'per') 
    eps = ob1_pad_eps(eps, t_pml * [1 1 0 0]);
end
% New size of padded epsilon structure
dims = size(eps);
pdims = prod(dims);
padded_height = dims(1);
padded_width = dims(2);
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
if strcmp(bc, 'pml')
    Ecurl = [scy(.5,.5)*-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
    Hcurl = [scy(.5,0)*(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 
elseif strcmp(bc, 'per')
    Ecurl = [-(S(0,1)-S(0,0)), scx(.5,.5)*(S(1,0)-S(0,0))];  
    Hcurl = [(S(0,0)-S(0,-1));  scx(0,.5)*-(S(0,0)-S(-1,0))]; 
end

% Diagonal matrix for 1/epsilon.
inv_eps = spdiags([eps_x(:).^-1; eps_y(:).^-1], 0, 2*prod(dims), 2*prod(dims));

% This is the matrix that we will solve.
A = Ecurl * inv_eps * Hcurl - omega^2 * speye(prod(dims));

%
% Determine the input excitation.
%

b = zeros(dims); % Input excitation, equivalent to magnetic current source.
in_pos = t_pml+2; % Cannot be 1, because eps interpolation wreaks havoc at border.

% For one-way excitation in the forward (to the right) direction,
% we simple cancel out excitation in the backward (left) direction.
if strcmp(bc, 'pml')
    b_pad = [floor((dims(2) - length(in.Hz))/2), ...
            ceil((dims(2) - length(in.Hz))/2)];
elseif strcmp(bc, 'per')
    b_pad = [0 0];
end
b(in_pos+1, b_pad(1)+1:end-b_pad(2)) = in.Hz;
b(in_pos, b_pad(1)+1:end-b_pad(2)) = -in.Hz * exp(1i * in.beta);

% Convert from field to current source.
b = b ./ eps_y;

% Vectorize.
b = b(:);

%
% Solve. This should be using sparse matrix factorization. 
%
Hz = A \ b;

%
% This piece is important for signifying how the objective fields are
% selected from the two-dimensional output field structure.  It is an input
% to the partial derivative with respect to the fields function handle
% because it ensures we put the adjoint source in the correct location.
%
S_filter = zeros(width, padded_width * padded_height);
compute_first = padded_width * padded_height - t_pml * (padded_height + 1);
for iter = 1 : 1 : width
    S_filter(width - (iter - 1), compute_first - (iter - 1) * padded_height) = 1.0;
end

vect_eps = eps(:);
complex_gradient = zeros(prod(dims), 1);

% Perform this multiplication once first instead of in every loop
Hcurl_Hz = Hcurl * Hz;

adj_src = g_x_partial(Hz, S_filter).';
% Vectorize to fit in form of the forward problem
adj_src = adj_src(:);
lambda = A' \ conj(adj_src);
lambda_adjoint = lambda';

% There are two parts to the total derivative.  The adjoint piece computed
% above as well as the partial derivative of the solution matrix A (in Ax =
% b) with respect to the parameter epsilon.
for e_col = 1 : 1 : width
    for e_row = 1 : 1 : height

        n = ((e_col - 1) + t_pml) * dims(1) + e_row + t_pml;

        get_eps_i = vect_eps(n);
        inv_eps2 = 0.5 / (get_eps_i^2);

        sparse_eps = zeros(2 * pdims, 1);

        sparse_eps(n) = inv_eps2 * Hcurl_Hz(n);
        sparse_eps(n + pdims) = inv_eps2 * Hcurl_Hz(n + pdims);

        % These two sets of if/else blocks take care of wrapping behavior
        % of epsilon interpolation in x and y.  For further optimization,
        % these blocks could be removed and the special cases at the
        % boundaries could be put in their own loops.
        if (mod(n, dims(1)) == 1)
           sparse_eps(n + dims(1) - 1) = inv_eps2 * Hcurl_Hz(n + dims(1) - 1);
        else
           sparse_eps(n - 1) = inv_eps2 * Hcurl_Hz(n - 1);
        end

        if (n <= dims(1))
            sparse_eps(n + pdims + (dims(2) - 1)*dims(1)) = inv_eps2 * Hcurl_Hz(n + pdims + (dims(2) - 1)*dims(1));
        else
            sparse_eps(n + pdims - dims(1)) = inv_eps2 * Hcurl_Hz(n + pdims - dims(1));
        end

        complex_gradient(n,1) = lambda_adjoint * (Ecurl * sparse_eps);

    end
end

% Obtain E-fields from Hz
E = (1i/omega) * inv_eps * Hcurl_Hz;

% Reshape and extract all three fields.
Ex = reshape(E(1:prod(dims)), dims);
Ey = reshape(E(prod(dims)+1:end), dims);
Hz = reshape(Hz, dims);

% Here, we take twice the real part because the variable we are taking
% derivative against is complex.
gradient = 2 * real(complex_gradient);
gradient = reshape(gradient, dims);

% Cut off the pml parts.
if strcmp(bc, 'pml')
    Ex = Ex(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
    Ey = Ey(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
    Hz = Hz(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
    gradient = gradient(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
elseif strcmp(bc, 'per')
    Ex = Ex(t_pml+1:end-t_pml, :);
    Ey = Ey(t_pml+1:end-t_pml, :);
    Hz = Hz(t_pml+1:end-t_pml, :);
    gradient = gradient(t_pml+1:end-t_pml, t_pml+1:end-t_pml);
end

