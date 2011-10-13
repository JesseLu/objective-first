function [eps] = solve(spec, max_iters, min_grad, varargin)
% EPS = SOLVE(SPEC, MAX_ITERS, MIN_GRAD, [METHOD])
%
% Description
%     Perform an objective-first optimization for a nanophotonic waveguide
%     coupler.
% 
% Inputs
%     SPEC: Structure.
%         The output of SETUP(). SPEC defines the waveguide coupler problem.
% 
%     MAX_ITERS: Non-negative integer.
%         Maximum number of iterations for which to run the optimization.
% 
%     MIN_GRAD: Non-negative scalar.
%         Minimum value for the norm of the gradient. The optimization will
%         stop once the norm of the gradient is below MIN_GRAD.
% 
%     METHOD: Character string (optional).
%         The name of the numerical method to use to solve the problem.
%         Currently only the 'sep_convex' option is supported.
% 
% Outputs
%     EPS: 2d array.
%         The permittivity values of the final design.
    
% Implementation notes
% *   We attempt to completely use linear-algebra language in this function,
%     which means everything is either a matrix or a vector (no 2-d arrays
%     for example). The implementation is much clearer this way.
% *   We fix the range of p to be from 0 to 1 (inclusive). In order to match
%     the desired range in epsilon values, p is scaled and offset.
% *   For numerical reasons (convexity in the p variable) we use the inverse
%     of epsilon instead of epsilon itself.


dims = size(spec.eps0);
N = prod(dims);


    % 
    % Form relevant matrices, and get initial values of x and p.
    %

[A, S] = ob1_matrices(dims); % Relevant matrices.

x0 = spec.Hz0(:); % Initial value of x (not needed for sep_convex method)
x_int = S.int * x0; % Interior values of the field (variable).
x_bnd = S.bnd * x0; % Boundary values of the field (fixed).

p_scale = diff(spec.eps_lims.^-1); % Scaling factor for p.
p_offset = (spec.eps_lims(1).^-1) * ones(N,1); % Offset for values of p.

p0 = p_scale^-1 * (spec.eps0(:).^-1 - p_offset); % Initial value of p.
p_int = S.int * p0; % Interior values of the structure (variable).
p_bnd = S.bnd * p0; % Boundary values of the structure (fixed).


    %
    % Construct helper functions.
    %

my_diag = @(z) spdiags(z, 0, length(z), length(z)); % Make diagonal matrix.

p2ie = @(p) p_scale * p + p_offset; % Transform from p to 1/eps.

% Physics operator in terms of x.
A_x = @(p) A{1} * my_diag(A{2} * p2ie(p)) * A{3} - spec.omega^2 * speye(N);

% Physics operator in terms of p (A_p * p - b_p = 0).
A_p = @(x) p_scale * A{1} * my_diag(A{3} * x) * A{2};
b_p = @(x) spec.omega^2 * x - A_p(x) * (p_offset ./ p_scale); 

% Physics residual.
phys_res = @(x, p) norm(A_x(p) * x)^2;
phys_res_alt = @(x, p) norm(A_p(x) * p - b_p(x))^2;

phys_res(x0, p0)
phys_res_alt(x0, p0)
eps = nan;
