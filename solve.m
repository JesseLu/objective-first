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
%         Note that the norm of the gradient is used, and not the norm-squared.
% 
%     METHOD: Character string (optional).
%         The name of the numerical method to use to solve the problem.
%         Currently only the 'alt_dir' option is supported. 'alt_dir'
%         requires the CVX package (www.stanford.edu/~boyd/cvx).
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

% Make sure we have access to the cvx package.
path(path, genpath('cvx'));

dims = size(spec.eps0);
N = prod(dims);


    % 
    % Form relevant matrices, and get initial values of x and p.
    %

[A, S] = ob1_matrices(dims); % Relevant matrices.

x0 = spec.Hz0(:); % Initial value of x.
x_int = S.int * x0; % Interior values of the field (variable).
x_bnd = S.bnd * x0; % Boundary values of the field (fixed).

p_scale = diff(spec.eps_lims.^-1); % Scaling factor for p.
p_offset = (spec.eps_lims(1).^-1) * ones(N,1); % Offset for values of p.

p0 = p_scale^-1 * (spec.eps0(:).^-1 - p_offset); % Initial value of p.
p_int = S.int * p0; % Interior values of the structure (variable).
p_bnd = S.bnd * p0; % Boundary values of the structure (fixed).


    %
    % Construct helper functions related to the physics residual.
    %

my_diag = @(z) spdiags(z, 0, length(z), length(z)); % Make diagonal matrix.

p2ie = @(p) p_scale * p + p_offset; % Transform from p to 1/eps.

% Physics operator in terms of x.
A_x = @(p) A{1} * my_diag(A{2} * p2ie(p)) * A{3} - spec.omega^2 * speye(N);

% Physics operator in terms of p (A_p * p - b_p = 0).
A_p = @(x) p_scale * A{1} * my_diag(A{3} * x) * A{2};
b_p = @(x) spec.omega^2 * x - A_p(x) * (p_offset ./ p_scale); 

% Physics residual.
phys_res = @(x, p) norm(S.res * (A_x(p) * x))^2;
% phys_res_alt = @(x, p) norm(S.res * (A_p(x) * p - b_p(x)))^2;

% Gradient of the physics residual. 
grad = @(x, p) [((S.res * A_x(p))' * (S.res * A_x(p)) * x); ...
                (S.res * A_p(x))'*(S.res * (A_p(x) * p - b_p(x)))];
                

    %
    % Helper functions for optimization process.
    %

% Functions to include add boundary values to interior values of x and p.
p_full = @(pin) S.int' * pin + S.bnd' * p_bnd;
x_full = @(xin) S.int' * xin + S.bnd' * x_bnd;

% Progress functions.
progress = @(xin, pin) [phys_res(x_full(xin), p_full(pin)), ...
                        norm(grad(x_full(xin), p_full(pin)))];
print_prog = @(iter, xin, pin) ...
    fprintf('%5d: \t%1.5e \t%1.5e\n', iter, progress(xin, pin));
print_header = @(str) ...
    fprintf('%s iter#: \t[phys_res] \t[grad_norm]\n', str);


    %
    % Perform the optimization.
    %

% Resolve what optimization method we should use.
if isempty(varargin)
    method = 'alt_dir';
else
    method = varargin{1};
end

switch method
    case 'alt_dir'
        % Alternating directions method.
        % *   The basic idea is to alternately optimize for x and then p.
        % *   This allows us to use standard, guaranteed-to-work solvers.
        % *   Although this method is slow, it should always approach a 
        %     solution.

        fprintf('Starting the alternating directions solver...\n');
        print_header('x/p'); % Print header information for progress.
        start_time = tic; % For timing purposes.

        for k = 1 : max_iters
            % Solve for x_int.
            x_int = (A_x(p_full(p_int)) * S.int') \ ...
                    (-A_x(p_full(p_int)) * S.bnd' * x_bnd);
            fprintf('(x) '); print_prog(k, x_int, p_int)

            % Solve for p_int.
            cvx_quiet(true);
            cvx_begin
                variable p_int(length(p_int))
                minimize norm(A_p(x_full(x_int)) * ...
                                (S.int' * p_int + S.bnd' * p_bnd) - ...
                                b_p(x_full(x_int)))
                subject to
                    p_int >= 0
                    p_int <= 1
            cvx_end
            fprintf('(p) '); print_prog(k, x_int, p_int)

            % Visualize.
            ob1_plot(dims,  {'p', p_full(p_int)}, ...
                            {'|x|', abs(x_full(x_int))}, ...
                            {'Re(x)', real(x_full(x_int))});

            % Check for gradient norm stopping condition.
            if (norm(grad(x_full(x_int), p_full(p_int))) < min_grad)
                fprintf('Gradient norm stopping condition satisfied.');
                break
            end
        end

    case 'alt_dir_mod'
        % Alternating directions method.
        % *   The basic idea is to alternately optimize for x and then p.
        % *   This allows us to use standard, guaranteed-to-work solvers.
        % *   Although this method is slow, it should always approach a 
        %     solution.

        fprintf('Starting the alternating directions solver...\n');
        print_header('x/p'); % Print header information for progress.
        start_time = tic; % For timing purposes.
        ptot = norm(p_int, 1);
        for k = 1 : max_iters
            % Solve for x_int.
            x_int = (A_x(p_full(p_int)) * S.int') \ ...
                    (-A_x(p_full(p_int)) * S.bnd' * x_bnd);
            fprintf('(x) '); print_prog(k, x_int, p_int)

            % Solve for p_int.
            cvx_quiet(true);
            cvx_begin
                variable p_int(length(p_int))
                minimize norm(A_p(x_full(x_int)) * ...
                                (S.int' * p_int + S.bnd' * p_bnd) - ...
                                b_p(x_full(x_int)))
                subject to
                    p_int >= 0
                    p_int <= 1
                    norm(p_int, 1) <= ptot
            cvx_end
            fprintf('(p) '); print_prog(k, x_int, p_int)

            % Visualize.
            ob1_plot(dims,  {'p', p_full(p_int)}, ...
                            {'|x|', abs(x_full(x_int))}, ...
                            {'Re(x)', real(x_full(x_int))});

            % Check for gradient norm stopping condition.
            if (norm(grad(x_full(x_int), p_full(p_int))) < min_grad)
                fprintf('Gradient norm stopping condition satisfied.');
                break
            end
        end

    otherwise
        error('Invalid choice of METHOD (%s).', method);
end

    
    %
    % Extract the full design and print final result.
    %

eps = reshape(p2ie(p_full(p_int)).^-1, dims); % Get full structure.

fprintf('\n%d iterations completed in %1.1f seconds.\n', k, toc(start_time));
fprintf('Final physics residual: %1.5e\nFinal gradient norm: %1.5e\n', ...
    progress(x_int, p_int));

% Visualize.
ob1_plot(dims,  {'p', p_full(p_int)}, ...
                {'|x|', abs(x_full(x_int))}, ...
                {'Re(x)', real(x_full(x_int))});

